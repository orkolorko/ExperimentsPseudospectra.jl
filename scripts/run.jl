using Distributed, ClusterManagers
procs = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))

@everywhere using Pkg; 
@everywhere Pkg.activate(@__DIR__)
@everywhere Pkg.instantiate(); 
@everywhere Pkg.precompile()


@info ENV["JULIA_CPU_TARGET"]
@everywhere ENV["OPENBLAS_NUM_THREADS"]=1

#procs = addprocs(2)

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

Ntot = ExperimentsPseudospectra.compute_steps(ρ, r_pearl)
@info "$N svd need to be computed"

start = 0
stop = 1000

@info "Start", start, "stop", stop

@async ExperimentsPseudospectra.submit_job(λ, ρ, r_pearl, jobs; start = start, stop = stop)

foreach(
      pid -> remote_do(ExperimentsPseudospectra.dowork, pid, D["P"], jobs, results),
      workers()[2:end]
)

#using DataFrames, JLD
#d = DataFrame()

avg_time = 0.0
N = stop-start

@elapsed while N > 0 # print out results
      x = take!(results)
      @info N, x

#      global avg_time += x.t
#      push!(d, x)

#     #  if N % 10000 == 0
#     #      @info "Done ", (1-Float64(N)/Ntot)*100, "%"
#     #  end 

      global N = N - 1
end
#save("results_$(λ)_$(ρ)_$(r_pearl)_$N_$(start)_$(stop).jld", "DataF", d)

avg_time /= stop-start

nworkers = length(procs)

@info "Average time for certifying an SVD", avg_time
# #@info "Estimate of time, serial" ceil(
# #    Int64, (2 * pi * 0.001) / (4.9 * 10^(-9)) * avg_time / 3600) " hours"
# #@info "Estimate of time, parallel, this configuration" ceil(
# #    Int64, (2 * pi * 0.001) / (4.9 * 10^(-9)) * avg_time / (3600 * nworkers)) " hours"

rmprocs(procs)