using Distributed, ClusterManagers, FileIO

procs = []

if haskey(ENV, "SLURM_NTASKS")
        procs = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
else
        procs = addprocs(4)
end


@everywhere using Pkg;
@everywhere Pkg.activate(@__DIR__)
@everywhere Pkg.instantiate();
@everywhere Pkg.precompile()

#@info ENV["JULIA_CPU_TARGET"]
@everywhere ENV["OPENBLAS_NUM_THREADS"] = 1

#
@everywhere using ExperimentsPseudospectra
@everywhere using JLD, CSV
@everywhere using LinearAlgebra
#@everywhere BLAS.set_num_threads(8)

@everywhere using BallArithmetic

@everywhere D = load("../ArnoldMatrixSchur256.jld")

const jobs = RemoteChannel(() -> Channel{Tuple}(32))
const results = RemoteChannel(() -> Channel{NamedTuple}(32))

S = D["S"]
ρ = 0.001
r_pearl = (4.9 * 10^(-5))
λ = 1.0

@info "Certifying ", λ, "radius", ρ, "radius pearl", r_pearl

start_angle = 0
stop_angle = 2 * pi

@everywhere include("script_functions.jl")

compute_enclosure_arc(
    D, λ, ρ, r_pearl; start_angle = start_angle,
    stop_angle = stop_angle, csvfile = "Arnold_one.csv")

rmprocs(procs)