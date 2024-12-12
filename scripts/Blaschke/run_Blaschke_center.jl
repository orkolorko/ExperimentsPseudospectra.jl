using Distributed, ClusterManagers, FileIO
using FileIO, Dates, CSV
#procs = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
procs = addprocs(2)

@everywhere using Pkg;
@everywhere Pkg.activate("../")
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

@everywhere D = load("./BlashkeMatrixSchur512.jld")

include("../script_functions.jl")

S = D["S"]
λ = 0.0


eig1 = S.values[1]
eig2 = S.values[2]

# we identify two sectors near the eigenvalues where we are going to work
# with a smaller step 

size_sec = 0.0001
angle_1 = angle(eig1)
angle_2 = 2 * pi + angle(eig2)

sectors = [(0.0, angle_1 - size_sec, 2^(-23))
           (angle_1 - size_sec, angle_1 + size_sec, 2^(-24));
           (angle_1 + size_sec, angle_2 - size_sec, 2^(-23));
           (angle_2 - size_sec, angle_2 + size_sec, 2^(-24));
           (angle_2 + size_sec, 2*pi, 2^(-23))]

for S in sectors
    start = S[1]
    stop = S[2]
    r_pearl = S[3]
    compute_enclosure_arc(
        D, λ, ρ, r_pearl; start_angle = π / 2, stop_angle = angle_1 - size_sec)
end

rmprocs(procs)