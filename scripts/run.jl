import Pkg;
Pkg.activate("./")

using Distributed, ClusterManagers
#addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))

procs = addprocs(2)

@everywhere using ExperimentsPseudospectra
@everywhere using JLD
@everywhere using LinearAlgebra
#@everywhere BLAS.set_num_threads(8)

@everywhere using BallArithmetic


@everywhere D = load("../ArnoldMatrixSchur256.jld")

const jobs = RemoteChannel(() -> Channel{Tuple}(32))
const results = RemoteChannel(() -> Channel{NamedTuple}(32))

ρ = 0.1
r_pearl = 0.01
λ = 1.0

N = ExperimentsPseudospectra.compute_steps(ρ, r_pearl)
@info "$N svd need to be computed"

@async ExperimentsPseudospectra.submit_job(1, ρ, r_pearl, jobs)

foreach(
    pid -> remote_do(ExperimentsPseudospectra.dowork, pid, D["P"], jobs, results),
    workers()
)

using DataFrames, JLD
d = DataFrame()

avg_time = 0.0
Ntot = N

@elapsed while N > 0 # print out results
    x = take!(results)

    global avg_time += x.t
    push!(d, x)

    global N = N - 1
end
save("results.jld", "DataF", d)

avg_time /= Ntot

nworkers = length(procs)

@info "Average time for certifying an SVD", avg_time
@info "Estimate of time, serial" ceil(
    Int64, (2 * pi * 0.001) / (4.9 * 10^(-9)) * avg_time / 3600) " hours"
@info "Estimate of time, parallel, this configuration" ceil(
    Int64, (2 * pi * 0.001) / (4.9 * 10^(-9)) * avg_time / (3600 * nworkers)) " hours"

rmprocs(procs)