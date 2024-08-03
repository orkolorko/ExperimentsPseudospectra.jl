import Pkg; 
Pkg.activate("./")

using Distributed, ClusterManagers
#addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))

procs = addprocs(2)

@everywhere using ExperimentsPseudospectra
@everywhere using JLD
@everywhere using BallArithmetic, LinearAlgebra

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

@elapsed while N > 0 # print out results
     x = take!(results)
     println(x)
     global N = N - 1
end



rmprocs(procs)