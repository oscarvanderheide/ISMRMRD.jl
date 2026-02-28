using HDF5
import HDF5.API: H5S_ALL, H5P_DEFAULT

# ---------------------------------------------------------------------------
# Memory layout constants for a 372-byte compound acquisition record
#
#   [0,   339] head  — AcquisitionHeader (340 bytes)
#   [340, 355] traj  — hvl_t{Float32}:  8-byte Csize_t len + 8-byte Ptr
#   [356, 371] data  — hvl_t{Float32}:  8-byte Csize_t len + 8-byte Ptr
# ---------------------------------------------------------------------------
const _COMPOUND_SIZE  = 372
const _TRAJ_HVL_OFF   = 340
const _DATA_HVL_OFF   = 356

# Absolute byte offsets of filter fields within the compound record
# (EncodingCounters starts at head offset 242)
const _OFF_SLICE       = 248  # 242 + 6
const _OFF_CONTRAST    = 250  # 242 + 8
const _OFF_REPETITION  = 254  # 242 + 12

# ---------------------------------------------------------------------------
# Public types
# ---------------------------------------------------------------------------

"""
    ISMRMRDFile(filename)

Reference to an ISMRMRD HDF5 file. Pass to `read_profiles`.
"""
struct ISMRMRDFile
    filename :: String
end

# ---------------------------------------------------------------------------
# VLen memory reclaim
# HDF5 1.12 renamed H5Dvlen_reclaim → H5Treclaim; handle both.
# ---------------------------------------------------------------------------
function _vlen_reclaim!(dtype_id, space_id, buf)
    if isdefined(HDF5.API, :h5t_reclaim)
        HDF5.API.h5t_reclaim(dtype_id, space_id, H5P_DEFAULT, buf)
    elseif isdefined(HDF5.API, :h5d_vlen_reclaim)
        HDF5.API.h5d_vlen_reclaim(dtype_id, space_id, H5P_DEFAULT, buf)
    else
        @warn "ISMRMRD: cannot locate VLen reclaim function in HDF5.API — memory may leak"
    end
end

# ---------------------------------------------------------------------------
# Filter helpers
# ---------------------------------------------------------------------------
_to_filter(::Nothing) = nothing
_to_filter(x::Integer) = (UInt16(x),)
_to_filter(xs) = Tuple(UInt16(x) for x in xs)

# Return a Ptr{UInt8} to the start of record i (0-based) in the flat buffer.
@inline function _rec(buf::Vector{UInt8}, i::Int)
    pointer(buf) + i * _COMPOUND_SIZE
end

# Build a Bool mask over N records using raw byte reads — no struct allocation.
function _filter_mask(buf::Vector{UInt8}, N::Int,
                      slice_f, contrast_f, repetition_f)
    mask = trues(N)
    GC.@preserve buf begin
        for i in 0:(N-1)
            p = _rec(buf, i)
            if !isnothing(slice_f)
                unsafe_load(Ptr{UInt16}(p + _OFF_SLICE)) ∉ slice_f &&
                    (@inbounds mask[i+1] = false; continue)
            end
            if !isnothing(contrast_f)
                unsafe_load(Ptr{UInt16}(p + _OFF_CONTRAST)) ∉ contrast_f &&
                    (@inbounds mask[i+1] = false; continue)
            end
            if !isnothing(repetition_f)
                unsafe_load(Ptr{UInt16}(p + _OFF_REPETITION)) ∉ repetition_f &&
                    (@inbounds mask[i+1] = false; continue)
            end
        end
    end
    return mask
end

# ---------------------------------------------------------------------------
# Profile construction from raw buffer
# ---------------------------------------------------------------------------

# Called from a GC.@preserve buf region; i is 0-based.
function _build_profile(buf::Vector{UInt8}, i::Int)
    p    = _rec(buf, i)
    head = parse_header(p)

    # Traj VLen: hvl_t at _TRAJ_HVL_OFF → (len::Csize_t, ptr::Ptr{Float32})
    traj_len = unsafe_load(Ptr{Csize_t}(p + _TRAJ_HVL_OFF))
    traj = if traj_len > 0
        traj_ptr = unsafe_load(Ptr{Ptr{Float32}}(p + _TRAJ_HVL_OFF + sizeof(Csize_t)))
        D   = Int(head.trajectory_dimensions)
        arr = unsafe_wrap(Vector{Float32}, traj_ptr, Int(traj_len))
        reshape(copy(arr), D, :)      # copy before vlen_reclaim frees traj_ptr
    else
        Matrix{Float32}(undef, 0, 0)
    end

    # Data VLen: hvl_t at _DATA_HVL_OFF → interleaved Re/Im for all channels
    data_len = unsafe_load(Ptr{Csize_t}(p + _DATA_HVL_OFF))
    data = if data_len > 0
        data_ptr = unsafe_load(Ptr{Ptr{Float32}}(p + _DATA_HVL_OFF + sizeof(Csize_t)))
        nsamples = Int(head.number_of_samples)
        nchans   = Int(head.active_channels)
        arr = unsafe_wrap(Vector{Float32}, data_ptr, Int(data_len))
        reshape(reinterpret(ComplexF32, copy(arr)), nsamples, nchans)
    else
        Matrix{ComplexF32}(undef, 0, 0)
    end

    return Profile(head, traj, data)
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    read_profiles(f::ISMRMRDFile;
                  dataset    = "dataset",
                  slice      = nothing,
                  contrast   = nothing,
                  repetition = nothing)
    -> RawAcquisitionData

Read all acquisitions from an ISMRMRD file using a single bulk HDF5 read
(O(1) HDF5 API calls regardless of the number of acquisitions).

Filtering by `slice`, `contrast`, or `repetition` is applied in Julia after
the read. Each filter accepts an integer, any iterable of integers, or
`nothing` (= keep all).

Returns a `RawAcquisitionData` with fields:
- `params`  — `Dict{String,Any}` parsed from the ISMRMRD XML header (same
               keys as MRIReco.jl's `GeneralParameters()`)
- `profiles` — `Vector{Profile}`

Profile construction is multithreaded; start Julia with multiple threads
(`julia -t auto`) to exploit this.
"""
function read_profiles(f::ISMRMRDFile;
                       dataset    = "dataset",
                       slice      = nothing,
                       contrast   = nothing,
                       repetition = nothing)

    slice_f      = _to_filter(slice)
    contrast_f   = _to_filter(contrast)
    repetition_f = _to_filter(repetition)

    h5open(f.filename, "r") do h
        xml_str = read(h["/$(dataset)/xml"])[1]
        dset    = h["/$(dataset)/data"]
        N       = prod(size(dset))

        # Flat byte buffer: N records × 372 bytes each.
        # HDF5 will fill VLen pointers (at offsets 340, 356 of each record)
        # with C-heap allocations that must later be freed via vlen_reclaim.
        buf   = Vector{UInt8}(undef, N * _COMPOUND_SIZE)
        dtype = _make_acquisition_type()

        # ── Single bulk read ────────────────────────────────────────────────
        HDF5.API.h5d_read(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf)

        # ── Vectorized filter on raw bytes (no struct allocation yet) ───────
        mask     = _filter_mask(buf, N, slice_f, contrast_f, repetition_f)
        idx_keep = findall(mask)
        profiles = Vector{Profile}(undef, length(idx_keep))

        # ── Threaded profile construction ───────────────────────────────────
        # buf stays live for the duration of @threads (outer task is blocked).
        GC.@preserve buf begin
            Threads.@threads for k in eachindex(idx_keep)
                profiles[k] = _build_profile(buf, idx_keep[k] - 1)
            end
        end

        # ── Free HDF5-allocated VLen heap memory ────────────────────────────
        space = HDF5.dataspace(dset)
        _vlen_reclaim!(dtype, space.id, buf)
        HDF5.API.h5t_close(dtype)
        close(space)

        return RawAcquisitionData(parse_params(xml_str), profiles)
    end
end
