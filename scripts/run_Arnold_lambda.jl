using Distributed, ClusterManagers, FileIO
#procs = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
procs = addprocs(2)

@everywhere using Pkg;
@everywhere Pkg.activate(@__DIR__)
@everywhere Pkg.instantiate();
@everywhere Pkg.precompile()

#@info ENV["JULIA_CPU_TARGET"]
@everywhere ENV["OPENBLAS_NUM_THREADS"] = 1

#
@everywhere using ExperimentsPseudospectra
@everywhere using JLD
@everywhere using LinearAlgebra
#@everywhere BLAS.set_num_threads(8)

@everywhere using BallArithmetic

@everywhere D = load("./ArnoldMatrixSchur256.jld")

const jobs = RemoteChannel(() -> Channel{Tuple}(32))
const results = RemoteChannel(() -> Channel{NamedTuple}(32))

S = D["S"]
ρ = 0.001
r_pearl = (4.9 * 10^(-9))
λ = S.values[1]

@info "Certifying ", λ, "radius", ρ, "radius pearl", r_pearl

start_angle = 0
stop_angle = 2*pi

compute_enclosure_arc(
        D, λ, ρ, r_pearl; start_angle = start_angle, stop_angle = stop_angle, csvfile = "Arnold_lambda.csv")

rmprocs(procs)