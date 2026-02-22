"""
    EncodingCounters

Encoding loop counters for one ISMRMRD acquisition. All fields are 0-indexed per the spec.
Embedded at byte offset 242 within AcquisitionHeader.
"""
struct EncodingCounters
    kspace_encode_step_1 :: UInt16
    kspace_encode_step_2 :: UInt16
    average              :: UInt16
    slice                :: UInt16
    contrast             :: UInt16
    phase                :: UInt16
    repetition           :: UInt16
    set                  :: UInt16
    segment              :: UInt16
    user                 :: NTuple{8, UInt16}
end

"""
    AcquisitionHeader

340-byte packed header describing one ISMRMRD acquisition (readout).
Field byte offsets are fixed by the ISMRMRD spec; see CLAUDE.md for the full table.

Note: `flags` is a UInt64 at byte offset 2 — misaligned relative to natural UInt64
alignment. Never attempt to reinterpret the whole struct from raw bytes.
"""
struct AcquisitionHeader
    version                :: UInt16
    flags                  :: UInt64
    measurement_uid        :: UInt32
    scan_counter           :: UInt32
    acquisition_time_stamp :: UInt32
    physiology_time_stamp  :: NTuple{3, UInt32}
    number_of_samples      :: UInt16
    available_channels     :: UInt16
    active_channels        :: UInt16
    channel_mask           :: NTuple{16, UInt64}
    discard_pre            :: UInt16
    discard_post           :: UInt16
    center_sample          :: UInt16
    encoding_space_ref     :: UInt16
    trajectory_dimensions  :: UInt16
    sample_time_us         :: Float32
    position               :: NTuple{3, Float32}
    read_dir               :: NTuple{3, Float32}
    phase_dir              :: NTuple{3, Float32}
    slice_dir              :: NTuple{3, Float32}
    patient_table_position :: NTuple{3, Float32}
    idx                    :: EncodingCounters
    user_int               :: NTuple{8, Int32}
    user_float             :: NTuple{8, Float32}
end

"""
    Profile

One acquired readout: header + optional k-space trajectory + complex data.

  - `traj`: (trajectory_dimensions × number_of_samples) or (0×0) if absent
  - `data`: (number_of_samples × active_channels)
"""
struct Profile
    head :: AcquisitionHeader
    traj :: Matrix{Float32}
    data :: Matrix{ComplexF32}
end

# ---------------------------------------------------------------------------
# Low-level byte-offset accessors
#
# All operate on a Ptr{UInt8} pointing to the start of the relevant field.
# Unaligned loads are hardware-safe on x86-64 and arm64 (Apple Silicon).
# ---------------------------------------------------------------------------

@inline rd_u16(p::Ptr{UInt8}, off::Int) = unsafe_load(Ptr{UInt16}(p + off))
@inline rd_u32(p::Ptr{UInt8}, off::Int) = unsafe_load(Ptr{UInt32}(p + off))
@inline rd_u64(p::Ptr{UInt8}, off::Int) = unsafe_load(Ptr{UInt64}(p + off))
@inline rd_f32(p::Ptr{UInt8}, off::Int) = unsafe_load(Ptr{Float32}(p + off))
@inline rd_i32(p::Ptr{UInt8}, off::Int) = unsafe_load(Ptr{Int32}(p + off))

@inline rd_ntuple_u16(p, off, ::Val{N}) where N = ntuple(i -> rd_u16(p, off + (i-1)*2), Val(N))
@inline rd_ntuple_u32(p, off, ::Val{N}) where N = ntuple(i -> rd_u32(p, off + (i-1)*4), Val(N))
@inline rd_ntuple_u64(p, off, ::Val{N}) where N = ntuple(i -> rd_u64(p, off + (i-1)*8), Val(N))
@inline rd_ntuple_f32(p, off, ::Val{N}) where N = ntuple(i -> rd_f32(p, off + (i-1)*4), Val(N))
@inline rd_ntuple_i32(p, off, ::Val{N}) where N = ntuple(i -> rd_i32(p, off + (i-1)*4), Val(N))

# p points to the start of the idx field (absolute offset 242 in the header)
function _parse_encoding_counters(p::Ptr{UInt8})
    EncodingCounters(
        rd_u16(p,  0),                  # kspace_encode_step_1
        rd_u16(p,  2),                  # kspace_encode_step_2
        rd_u16(p,  4),                  # average
        rd_u16(p,  6),                  # slice
        rd_u16(p,  8),                  # contrast
        rd_u16(p, 10),                  # phase
        rd_u16(p, 12),                  # repetition
        rd_u16(p, 14),                  # set
        rd_u16(p, 16),                  # segment
        rd_ntuple_u16(p, 18, Val(8)),   # user[8]
    )
end

"""
    parse_header(p::Ptr{UInt8}) -> AcquisitionHeader

Parse a 340-byte packed AcquisitionHeader from a raw byte pointer using
explicit byte offsets. Zero allocation, zero reflection.
"""
function parse_header(p::Ptr{UInt8})
    AcquisitionHeader(
        rd_u16(p,   0),                  # version
        rd_u64(p,   2),                  # flags  (misaligned — intentional)
        rd_u32(p,  10),                  # measurement_uid
        rd_u32(p,  14),                  # scan_counter
        rd_u32(p,  18),                  # acquisition_time_stamp
        rd_ntuple_u32(p, 22, Val(3)),    # physiology_time_stamp[3]
        rd_u16(p,  34),                  # number_of_samples
        rd_u16(p,  36),                  # available_channels
        rd_u16(p,  38),                  # active_channels
        rd_ntuple_u64(p, 40, Val(16)),   # channel_mask[16]
        rd_u16(p, 168),                  # discard_pre
        rd_u16(p, 170),                  # discard_post
        rd_u16(p, 172),                  # center_sample
        rd_u16(p, 174),                  # encoding_space_ref
        rd_u16(p, 176),                  # trajectory_dimensions
        rd_f32(p, 178),                  # sample_time_us
        rd_ntuple_f32(p, 182, Val(3)),   # position[3]
        rd_ntuple_f32(p, 194, Val(3)),   # read_dir[3]
        rd_ntuple_f32(p, 206, Val(3)),   # phase_dir[3]
        rd_ntuple_f32(p, 218, Val(3)),   # slice_dir[3]
        rd_ntuple_f32(p, 230, Val(3)),   # patient_table_position[3]
        _parse_encoding_counters(p + 242),  # idx (34 bytes)
        rd_ntuple_i32(p, 276, Val(8)),   # user_int[8]
        rd_ntuple_f32(p, 308, Val(8)),   # user_float[8]
    )
end
