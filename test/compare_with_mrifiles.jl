"""
Compare ISMRMRD.jl output against MRIReco/MRIFiles as a reference.

Usage:
    julia --project=/path/to/ISMRMRD test/compare_with_mrifiles.jl /path/to/file.h5

MRIFiles and MRIBase must be installed in the active environment (or the
default environment if they are not in the ISMRMRD project):

    julia -e 'import Pkg; Pkg.add(["MRIFiles","MRIBase"])'

What is checked:
  1. params dict — all keys present, values structurally equal
  2. profile count
  3. Per-profile: every AcquisitionHeader field, traj matrix, data matrix
"""

using Test

# ── Imports ───────────────────────────────────────────────────────────────────
# Use explicit module prefixes throughout to avoid name collisions between
# ISMRMRD.ISMRMRDFile and MRIFiles.ISMRMRDFile, etc.

import MRIFiles
import MRIBase
import ISMRMRD

# ── Path argument ─────────────────────────────────────────────────────────────

isempty(ARGS) && error("""
No file path supplied.
Usage: julia --project test/compare_with_mrifiles.jl /path/to/file.h5
""")

path = ARGS[1]
isfile(path) || error("File not found: $path")

# ── Read with both implementations ────────────────────────────────────────────

println("Reading with MRIFiles … ")
ref_raw = MRIBase.RawAcquisitionData(MRIFiles.ISMRMRDFile(path))

println("Reading with ISMRMRD.jl … ")
got_raw = ISMRMRD.read_profiles(ISMRMRD.ISMRMRDFile(path))

println()

# ── Structural equality helper ────────────────────────────────────────────────
# Compares two values that may come from identically-structured types in
# different modules (e.g. MRIBase.Limit vs ISMRMRD.Limit).
# Falls back to field-by-field comparison when types differ.

function structural_equal(a, b)
    # Fast path: same type, delegate to ==
    typeof(a) === typeof(b) && return a == b

    # Sequences: compare element-by-element
    if (a isa Union{Tuple, AbstractArray}) && (b isa Union{Tuple, AbstractArray})
        length(a) == length(b) || return false
        return all(structural_equal(x, y) for (x, y) in zip(a, b))
    end

    # Composite types: compare field-by-field
    fa = fieldnames(typeof(a))
    fb = fieldnames(typeof(b))
    isempty(fa) && return false   # not a composite type, already ≠
    fa == fb     || return false  # different structure
    return all(structural_equal(getfield(a, f), getfield(b, f)) for f in fa)
end

# ── 1. params ─────────────────────────────────────────────────────────────────

@testset "params" begin
    ref_keys = Set(keys(ref_raw.params))
    got_keys = Set(keys(got_raw.params))

    missing_in_got = setdiff(ref_keys, got_keys)
    extra_in_got   = setdiff(got_keys, ref_keys)

    if !isempty(missing_in_got)
        @warn "Keys present in MRIFiles but missing in ISMRMRD.jl" missing_in_got
    end
    if !isempty(extra_in_got)
        @warn "Keys present in ISMRMRD.jl but absent in MRIFiles" extra_in_got
    end

    @test isempty(missing_in_got)
    @test isempty(extra_in_got)

    for k in sort(collect(intersect(ref_keys, got_keys)))
        rv = ref_raw.params[k]
        gv = got_raw.params[k]
        ok = structural_equal(rv, gv)
        if !ok
            @warn "Mismatch" key=k ref=rv got=gv
        end
        @test ok broken=false
    end
end

# ── 2. profiles ───────────────────────────────────────────────────────────────

@testset "profile count" begin
    @test length(ref_raw.profiles) == length(got_raw.profiles)
end

n = min(length(ref_raw.profiles), length(got_raw.profiles))

@testset "profiles" begin
    @testset "profile $i" for i in 1:n
        rp = ref_raw.profiles[i]
        gp = got_raw.profiles[i]

        # AcquisitionHeader — compare every field by name
        @testset "head" begin
            for f in fieldnames(typeof(rp.head))
                rv = getfield(rp.head, f)
                gv = getfield(gp.head, f)
                ok = structural_equal(rv, gv)
                if !ok
                    @warn "Header field mismatch" profile=i field=f ref=rv got=gv
                end
                @test ok
            end
        end

        @testset "traj" begin
            @test size(rp.traj) == size(gp.traj)
            @test rp.traj == gp.traj
        end

        @testset "data" begin
            @test size(rp.data) == size(gp.data)
            @test rp.data == gp.data
        end
    end
end

println("\nDone.")
