# HDF5 compound type construction for ISMRMRD acquisitions.
# Mirrors the on-disk layout defined by the ISMRMRD spec.
# Adapted from MRIReco.jl/MRIFiles/src/ISMRMRD/HDF5Types.jl.

import HDF5.API: h5t_create, h5t_insert, h5t_array_create, h5t_close,
                 h5t_vlen_create, hsize_t, H5T_COMPOUND
import HDF5: hdf5_type_id

# Build the compound type for EncodingCounters (34 bytes).
function _make_encoding_counters_type()
    dt = h5t_create(H5T_COMPOUND, 34)
    h5t_insert(dt, "kspace_encode_step_1",  0, hdf5_type_id(UInt16))
    h5t_insert(dt, "kspace_encode_step_2",  2, hdf5_type_id(UInt16))
    h5t_insert(dt, "average",               4, hdf5_type_id(UInt16))
    h5t_insert(dt, "slice",                 6, hdf5_type_id(UInt16))
    h5t_insert(dt, "contrast",              8, hdf5_type_id(UInt16))
    h5t_insert(dt, "phase",                10, hdf5_type_id(UInt16))
    h5t_insert(dt, "repetition",           12, hdf5_type_id(UInt16))
    h5t_insert(dt, "set",                  14, hdf5_type_id(UInt16))
    h5t_insert(dt, "segment",              16, hdf5_type_id(UInt16))
    d = h5t_array_create(hdf5_type_id(UInt16), Cuint(1), hsize_t[8])
    h5t_insert(dt, "user", 18, d)
    h5t_close(d)
    return dt
end

# Build the compound type for AcquisitionHeader (340 bytes).
function _make_acquisition_header_type()
    dt = h5t_create(H5T_COMPOUND, 340)
    h5t_insert(dt, "version",                  0, hdf5_type_id(UInt16))
    h5t_insert(dt, "flags",                    2, hdf5_type_id(UInt64))
    h5t_insert(dt, "measurement_uid",         10, hdf5_type_id(UInt32))
    h5t_insert(dt, "scan_counter",            14, hdf5_type_id(UInt32))
    h5t_insert(dt, "acquisition_time_stamp",  18, hdf5_type_id(UInt32))

    d = h5t_array_create(hdf5_type_id(UInt32), Cuint(1), hsize_t[3])
    h5t_insert(dt, "physiology_time_stamp", 22, d); h5t_close(d)

    h5t_insert(dt, "number_of_samples",  34, hdf5_type_id(UInt16))
    h5t_insert(dt, "available_channels", 36, hdf5_type_id(UInt16))
    h5t_insert(dt, "active_channels",    38, hdf5_type_id(UInt16))

    d = h5t_array_create(hdf5_type_id(UInt64), Cuint(1), hsize_t[16])
    h5t_insert(dt, "channel_mask", 40, d); h5t_close(d)

    h5t_insert(dt, "discard_pre",           168, hdf5_type_id(UInt16))
    h5t_insert(dt, "discard_post",          170, hdf5_type_id(UInt16))
    h5t_insert(dt, "center_sample",         172, hdf5_type_id(UInt16))
    h5t_insert(dt, "encoding_space_ref",    174, hdf5_type_id(UInt16))
    h5t_insert(dt, "trajectory_dimensions", 176, hdf5_type_id(UInt16))
    h5t_insert(dt, "sample_time_us",        178, hdf5_type_id(Float32))

    d = h5t_array_create(hdf5_type_id(Float32), Cuint(1), hsize_t[3])
    h5t_insert(dt, "position",               182, d)
    h5t_insert(dt, "read_dir",               194, d)
    h5t_insert(dt, "phase_dir",              206, d)
    h5t_insert(dt, "slice_dir",              218, d)
    h5t_insert(dt, "patient_table_position", 230, d)
    h5t_close(d)

    denc = _make_encoding_counters_type()
    h5t_insert(dt, "idx", 242, denc); h5t_close(denc)

    d = h5t_array_create(hdf5_type_id(Int32),   Cuint(1), hsize_t[8])
    h5t_insert(dt, "user_int", 276, d); h5t_close(d)

    d = h5t_array_create(hdf5_type_id(Float32), Cuint(1), hsize_t[8])
    h5t_insert(dt, "user_float", 308, d); h5t_close(d)

    return dt
end

"""
    _make_acquisition_type() -> hid_t

Build the full 372-byte HDF5 compound type for one acquisition record:

  offset   0 : head (AcquisitionHeader, 340 bytes)
  offset 340 : traj (hvl_t{Float32}, 16 bytes — len::size_t + ptr::void*)
  offset 356 : data (hvl_t{Float32}, 16 bytes)

Caller is responsible for closing the returned type ID with `h5t_close`.
"""
function _make_acquisition_type()
    dt = h5t_create(H5T_COMPOUND, 372)

    dhead = _make_acquisition_header_type()
    h5t_insert(dt, "head", 0, dhead); h5t_close(dhead)

    # Two separate vlen types (each must be independently created and closed)
    dvlen = h5t_vlen_create(hdf5_type_id(Float32))
    h5t_insert(dt, "traj", 340, dvlen); h5t_close(dvlen)

    dvlen = h5t_vlen_create(hdf5_type_id(Float32))
    h5t_insert(dt, "data", 356, dvlen); h5t_close(dvlen)

    return dt
end
