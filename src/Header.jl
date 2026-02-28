using LightXML

# ---------------------------------------------------------------------------
# Supporting types  (mirror MRIBase definitions so downstream code is compatible)
# ---------------------------------------------------------------------------

struct Limit
    minimum :: Int
    maximum :: Int
    center  :: Int
end

struct MeasurementDependency
    dependencyType :: String
    measurementID  :: String
end

struct CoilDescription
    number :: Int
    name   :: String
end

"""
    RawAcquisitionData

Container returned by `read_profiles`, holding the parsed XML parameter
dictionary and the vector of acquired profiles.
Field names and `params` keys match those produced by MRIReco.jl so that
downstream code using `RawAcquisitionData` from MRIReco can be used unchanged.
"""
struct RawAcquisitionData
    params   :: Dict{String, Any}
    profiles :: Vector{Profile}
end

# ---------------------------------------------------------------------------
# Low-level XML helpers
# ---------------------------------------------------------------------------

# Add a value (or vector of values) to params under `key`, parsed as type T.
# Silent no-op when the tag is absent.
function _add_to_dict!(params, el, tag, ::Type{T}, key=tag) where T
    nodes = get_elements_by_tagname(el, tag)
    isempty(nodes) && return
    if length(nodes) == 1
        params[key] = _parse_node(T, nodes[1])
    else
        params[key] = [_parse_node(T, n) for n in nodes]
    end
end

_parse_node(::Type{T}, node) where T = parse(T, content(node))
_parse_node(::Type{String}, node)    = content(node)

function _parse_node(::Type{Vector{T}}, node) where T
    x = content(get_elements_by_tagname(node, "x")[1])
    y = content(get_elements_by_tagname(node, "y")[1])
    z = content(get_elements_by_tagname(node, "z")[1])
    return [parse(T, x), parse(T, y), parse(T, z)]
end

function _parse_node(::Type{Limit}, node)
    Limit(
        parse(Int, content(get_elements_by_tagname(node, "minimum")[1])),
        parse(Int, content(get_elements_by_tagname(node, "maximum")[1])),
        parse(Int, content(get_elements_by_tagname(node, "center")[1])),
    )
end

function _parse_node(::Type{MeasurementDependency}, node)
    MeasurementDependency(
        content(get_elements_by_tagname(node, "dependencyType")[1]),
        content(get_elements_by_tagname(node, "measurementID")[1]),
    )
end

function _parse_node(::Type{CoilDescription}, node)
    CoilDescription(
        parse(Int, content(get_elements_by_tagname(node, "coilNumber")[1])),
        content(get_elements_by_tagname(node, "coilName")[1]),
    )
end

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    parse_params(xml_str::String) -> Dict{String, Any}

Parse an ISMRMRD XML header string into a parameter dictionary whose keys and
structure match those produced by MRIReco.jl's `GeneralParameters()`.
"""
function parse_params(xml_str::String)
    xdoc = parse_string(xml_str)
    params = _extract_params(xdoc)
    free(xdoc)
    return params
end

# ---------------------------------------------------------------------------
# Internal XML extraction  (ported from MRIReco MRIFiles/src/ISMRMRD/Parameters.jl)
# ---------------------------------------------------------------------------

function _extract_params(xdoc::XMLDocument)
    params = Dict{String, Any}()
    root   = LightXML.root(xdoc)

    # SubjectInformation
    e = get_elements_by_tagname(root, "subjectInformation")
    if !isempty(e)
        _add_to_dict!(params, e[1], "patientName",      String)
        _add_to_dict!(params, e[1], "patientWeight_kg", Float64)
        _add_to_dict!(params, e[1], "patientID",        String)
        _add_to_dict!(params, e[1], "patientBirthdate", String)
        _add_to_dict!(params, e[1], "patientGender",    String)
    end

    # StudyInformation
    e = get_elements_by_tagname(root, "studyInformation")
    if !isempty(e)
        _add_to_dict!(params, e[1], "studyDate",              String)
        _add_to_dict!(params, e[1], "studyID",                String)
        _add_to_dict!(params, e[1], "accessionNumber",        Int)
        _add_to_dict!(params, e[1], "referringPhysicianName", String)
        _add_to_dict!(params, e[1], "studyDescription",       String)
        _add_to_dict!(params, e[1], "studyInstanceUID",       String)
    end

    # MeasurementInformation
    e = get_elements_by_tagname(root, "measurementInformation")
    if !isempty(e)
        _add_to_dict!(params, e[1], "measurementID",           String)
        _add_to_dict!(params, e[1], "seriesDate",              String)
        _add_to_dict!(params, e[1], "seriesTime",              String)
        _add_to_dict!(params, e[1], "patientPosition",         String)
        _add_to_dict!(params, e[1], "initialSeriesNumber",     Int)
        _add_to_dict!(params, e[1], "protocolName",            String)
        _add_to_dict!(params, e[1], "seriesDescription",       String)
        _add_to_dict!(params, e[1], "seriesInstanceUIDRoot",   String)
        _add_to_dict!(params, e[1], "frameOfReferenceUID",     String)
        _add_to_dict!(params, e[1], "measurementDependency",   MeasurementDependency)
        _add_to_dict!(params, e[1], "referencedImageSequence", String)
    end

    # AcquisitionSystemInformation
    e = get_elements_by_tagname(root, "acquisitionSystemInformation")
    if !isempty(e)
        _add_to_dict!(params, e[1], "systemVendor",                   String)
        _add_to_dict!(params, e[1], "systemModel",                    String)
        _add_to_dict!(params, e[1], "systemFieldStrength_T",          Float64)
        _add_to_dict!(params, e[1], "relativeReceiverNoiseBandwidth", Float64)
        _add_to_dict!(params, e[1], "receiverChannels",               Int)
        _add_to_dict!(params, e[1], "institutionName",                String)
        _add_to_dict!(params, e[1], "stationName",                    String)
        _add_to_dict!(params, e[1], "coilLabel",                      CoilDescription)
    end

    # ExperimentalConditions
    e = get_elements_by_tagname(root, "experimentalConditions")
    if !isempty(e)
        _add_to_dict!(params, e[1], "H1resonanceFrequency_Hz", Int)
    end

    # Encoding
    e = get_elements_by_tagname(root, "encoding")
    if !isempty(e)
        enc = e[1]

        es = enc["encodedSpace"]
        if !isempty(es)
            params["encodedSize"] = _parse_node(Vector{Int},     es[1]["matrixSize"][1])
            params["encodedFOV"]  = _parse_node(Vector{Float64}, es[1]["fieldOfView_mm"][1])
        end

        rs = enc["reconSpace"]
        if !isempty(rs)
            params["reconSize"] = _parse_node(Vector{Int},     rs[1]["matrixSize"][1])
            params["reconFOV"]  = _parse_node(Vector{Float64}, rs[1]["fieldOfView_mm"][1])
        end

        d = enc["encodingLimits"]
        if !isempty(d)
            _add_to_dict!(params, d[1], "kspace_encoding_step_0", Limit, "enc_lim_kspace_encoding_step_0")
            _add_to_dict!(params, d[1], "kspace_encoding_step_1", Limit, "enc_lim_kspace_encoding_step_1")
            _add_to_dict!(params, d[1], "kspace_encoding_step_2", Limit, "enc_lim_kspace_encoding_step_2")
            _add_to_dict!(params, d[1], "average",    Limit, "enc_lim_average")
            _add_to_dict!(params, d[1], "slice",      Limit, "enc_lim_slice")
            _add_to_dict!(params, d[1], "contrast",   Limit, "enc_lim_contrast")
            _add_to_dict!(params, d[1], "phase",      Limit, "enc_lim_phase")
            _add_to_dict!(params, d[1], "repetition", Limit, "enc_lim_repetition")
            _add_to_dict!(params, d[1], "set",        Limit, "enc_lim_set")
            _add_to_dict!(params, d[1], "segment",    Limit, "enc_lim_segment")
        end

        _add_to_dict!(params, enc, "trajectory", String)

        d = enc["trajectoryDescription"]
        if !isempty(d)
            y = Dict{String, Any}()
            y["identifier"] = content(d[1]["identifier"][1])
            for q in d[1]["userParameterLong"]
                y[content(q["name"][1])] = parse(Int,     content(q["value"][1]))
            end
            for q in d[1]["userParameterDouble"]
                y[content(q["name"][1])] = parse(Float64, content(q["value"][1]))
            end
            for q in d[1]["userParameterString"]
                y[content(q["name"][1])] = content(q["value"][1])
            end
            if !isempty(d[1]["comment"])
                y["comment"] = content(d[1]["comment"][1])
            end
            params["trajectoryDescription"] = y
        end

        d = enc["parallelImaging"]
        if !isempty(d)
            a = parse(Int, content(d[1]["accelerationFactor"][1]["kspace_encoding_step_1"][1]))
            b = parse(Int, content(d[1]["accelerationFactor"][1]["kspace_encoding_step_2"][1]))
            params["accelerationFactor"] = (a, b)
            _add_to_dict!(params, d[1], "calibrationMode",       String)
            _add_to_dict!(params, d[1], "interleavingDimension", String)
        end

        _add_to_dict!(params, enc, "echoTrainLength", Int)
    end

    # Support incorrectly-created files with bare <parallelImaging> at root level
    e = get_elements_by_tagname(root, "parallelImaging")
    if !isempty(e)
        a = parse(Int, content(e[1]["accelerationFactor"][1]["kspace_encoding_step_1"][1]))
        b = parse(Int, content(e[1]["accelerationFactor"][1]["kspace_encoding_step_2"][1]))
        params["accelerationFactor"] = (a, b)
        _add_to_dict!(params, e[1], "calibrationMode",       String)
        _add_to_dict!(params, e[1], "interleavingDimension", String)
    end

    # SequenceParameters
    e = get_elements_by_tagname(root, "sequenceParameters")
    if !isempty(e)
        _add_to_dict!(params, e[1], "TR",            Float64)
        _add_to_dict!(params, e[1], "TE",            Float64)
        _add_to_dict!(params, e[1], "TI",            Float64)
        _add_to_dict!(params, e[1], "flipAngle_deg", Float64)
        _add_to_dict!(params, e[1], "sequence_type", String)
        _add_to_dict!(params, e[1], "echo_spacing",  Float64)
    end

    # WaveformInformation
    e = get_elements_by_tagname(root, "waveformInformation")
    if !isempty(e)
        _add_to_dict!(params, e[1], "waveformName", Float64)
        _add_to_dict!(params, e[1], "waveformType", Float64)
        d = e[1]["userParameters"]
        if !isempty(d)
            y = Dict{String, Any}()
            for q in d[1]["userParameterLong"]
                y[content(q["name"][1])] = parse(Int,     content(q["value"][1]))
            end
            for q in d[1]["userParameterDouble"]
                y[content(q["name"][1])] = parse(Float64, content(q["value"][1]))
            end
            for q in d[1]["userParameterString"]
                y[content(q["name"][1])] = content(q["value"][1])
            end
            params["waveformUserParameters"] = y
        end
    end

    # UserParameters
    e = get_elements_by_tagname(root, "userParameters")
    if !isempty(e)
        y = Dict{String, Any}()
        for q in e[1]["userParameterLong"]
            y[content(q["name"][1])] = parse(Int,     content(q["value"][1]))
        end
        for q in e[1]["userParameterDouble"]
            y[content(q["name"][1])] = parse(Float64, content(q["value"][1]))
        end
        for q in e[1]["userParameterString"]
            y[content(q["name"][1])] = content(q["value"][1])
        end
        params["userParameters"] = y
    end

    return params
end
