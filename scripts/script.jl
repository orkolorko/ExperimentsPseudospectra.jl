using Distributed, ClusterManagers, FileIO
using FileIO, Dates, CSV
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

include("script_functions.jl")

S = D["S"]
ρ = 0.001
r_pearl = (4.9 * 10^(-9))
λ = S.values[1]

compute_enclosure_arc(D, λ, ρ, r_pearl; start_angle = 0.0, stop_angle = π/1000)

rmprocs(procs)