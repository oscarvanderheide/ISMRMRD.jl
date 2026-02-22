# ISMRMRD.jl

Fast Julia reader for the [ISMRMRD](https://ismrmrd.readthedocs.io/en/latest/) MRI raw data format.

Compared to the reader in [MRIReco.jl](https://github.com/MagneticResonanceImaging/MRIReco.jl), this package reads the entire acquisition dataset in a **single HDF5 call** (O(1) instead of O(N)), deserialises headers via direct byte-offset access (no `IOBuffer`, no reflection), and constructs profiles in parallel with `Threads.@threads`.

---

## Usage

```julia
using ISMRMRD

f = ISMRMRDFile("my_scan.h5")

# Read all acquisitions
profiles, xml = read_profiles(f)

# Filter by index (0-based, matching ISMRMRD convention)
profiles, xml = read_profiles(f; slice=0, contrast=0:2, repetition=5)

# Access profile data
p = profiles[1]
p.head.number_of_samples   # UInt16
p.head.active_channels      # UInt16
p.head.idx.slice            # UInt16
p.data                      # Matrix{ComplexF32}: (samples × channels)
p.traj                      # Matrix{Float32}:   (dims × samples), or (0×0)

# The raw XML parameter string from the file header
println(xml)
```

`read_profiles` returns a `Vector{Profile}` and the raw XML header string.
Start Julia with multiple threads (`julia -t auto`) to speed up profile construction.

### Exported types

| Type | Description |
|---|---|
| `ISMRMRDFile` | Thin wrapper around a filename |
| `AcquisitionHeader` | 340-byte header struct (all ISMRMRD fields) |
| `EncodingCounters` | Encoding loop counters embedded in the header |
| `Profile` | One readout: header + traj + data |

---

## Developer Guide

### Architecture

```
src/
  ISMRMRD.jl           module entry point
  AcquisitionHeader.jl struct definitions + byte-offset header parser
  HDF5Types.jl         HDF5 compound type construction
  Reader.jl            bulk VLen read, filtering, profile assembly
```

### Why it is fast

#### 1. Single bulk HDF5 read

The ISMRMRD `/dataset/data` HDF5 dataset is a 1-D array of compound records. Each record is 372 bytes:

```
[0,   339]  head: AcquisitionHeader (340 bytes, packed)
[340, 355]  traj: hvl_t{Float32}   (8-byte len + 8-byte ptr, C heap)
[356, 371]  data: hvl_t{Float32}   (8-byte len + 8-byte ptr, C heap)
```

MRIReco.jl reads each of the N acquisitions individually (`h[...][m]` in a loop = N separate HDF5 C-library calls). This package reads all N records in one call:

```julia
buf   = Vector{UInt8}(undef, N * 372)
dtype = _make_acquisition_type()
HDF5.API.h5d_read(dset, dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf)
```

After the read, `buf` holds all headers plus `hvl_t` structs whose pointers reference C-heap-allocated VLen data. The data is copied into Julia arrays and then the C heap is freed:

```julia
HDF5.API.h5d_vlen_reclaim(dtype, space.id, H5P_DEFAULT, buf)
```

#### 2. Zero-allocation header parsing

`AcquisitionHeader` has a packed, non-naturally-aligned layout (e.g. `flags::UInt64` sits at byte offset 2). Julia's struct layout adds padding, so `sizeof(AcquisitionHeader) == 352`, not 340. A direct `reinterpret` of raw bytes would give wrong results.

Instead, `parse_header(p::Ptr{UInt8})` reads each field by its known byte offset using `unsafe_load`:

```julia
@inline rd_u64(p::Ptr{UInt8}, off::Int) = unsafe_load(Ptr{UInt64}(p + off))
# ...
rd_u64(p, 2)   # flags at offset 2 — intentionally unaligned
```

Unaligned loads are safe (and fast) on x86-64 and arm64 (Apple Silicon). No `IOBuffer`, no `setfield!` reflection loop, no allocation per header.

#### 3. Vectorised pre-filtering

`slice`, `contrast`, and `repetition` filters are evaluated directly on the raw byte buffer — reading the relevant `UInt16` at known offsets — before any Julia struct is allocated. Only the matching acquisitions go through `_build_profile`.

#### 4. Threaded profile construction

Profile construction (copy VLen data from C heap, reinterpret as `ComplexF32`, reshape) is embarrassingly parallel and runs under `Threads.@threads`.

### VLen memory ownership

After `h5d_read`, the `hvl_t` pointers in `buf` point to memory owned by the HDF5 C library. The pattern is:

```
h5d_read  →  copy VLen data into Julia arrays  →  h5d_vlen_reclaim
```

`_build_profile` always calls `copy` on the `unsafe_wrap`-ed arrays before returning, so by the time `vlen_reclaim` runs, all data is safely in Julia-owned memory.

### HDF5.jl API compatibility

| Function | Available in |
|---|---|
| `HDF5.hdf5_type_id` | HDF5.jl ≥ 0.16 |
| `HDF5.API.h5d_vlen_reclaim` | HDF5.jl 0.16–0.17 |
| `HDF5.API.h5t_reclaim` | HDF5.jl ≥ 0.17 (HDF5 C lib ≥ 1.12) |

`_vlen_reclaim!` checks for both at runtime and uses whichever is available.

### Thread safety note

All HDF5 I/O runs on a single thread (inside `h5open`). Only the Julia-side processing (header parsing, data copying, struct construction) is parallelised. HDF5.jl is not safe for concurrent reads on the same file handle.
