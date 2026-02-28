module ISMRMRD

using HDF5
using LightXML

include("AcquisitionHeader.jl")
include("Header.jl")
include("HDF5Types.jl")
include("Reader.jl")

export ISMRMRDFile, read_profiles
export AcquisitionHeader, EncodingCounters, Profile
export RawAcquisitionData, Limit, MeasurementDependency, CoilDescription
export parse_params

end # module ISMRMRD
