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
λ = 0.0
ρ = 0.27812055400531616

eig1 = S.values[1]
eig2 = S.values[2]

# we identify two sectors near the eigenvalues where we are going to work
# with a smaller step 

size_sec = 0.0001
angle_1 = angle(eig1)
angle_2 = 2 * pi + angle(eig2)

sectors = [(angle_1 - size_sec, angle_1 + size_sec, 6.36351520069046e-8);
           (angle_2 - size_sec, angle_2 + size_sec, 6.36351520069046e-8)]

# r_pearl = 6.36351520069046e-8
# compute_enclosure_arc(D, λ, ρ, r_pearl; start_angle = angle_1 - size_sec, stop_angle = angle_1 + size_sec)
# 
# compute_enclosure_arc(D, λ, ρ, r_pearl; start_angle = angle_2 - size_sec, stop_angle = angle_2 + size_sec)

r_pearl = 0.0000001
compute_enclosure_arc(
    D, λ, ρ, r_pearl; start_angle = π / 2, stop_angle = angle_1 - size_sec)

rmprocs(procs)