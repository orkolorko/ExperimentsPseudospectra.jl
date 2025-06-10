import Pkg;
Pkg.activate(@__DIR__)

using Logging, Dates, Distributed, LinearAlgebra, ClusterManagers, DataFrames, JLD2, CSV

datetime = Dates.now()

if haskey(ENV, "SLURM_NTASKS")
    procs = addprocs_slurm(parse(Int, ENV["SLURM_NTASKS"]))
    location = "slurm"
else
    procs = addprocs(4)
    location = "local"
end

nprocs = length(procs)

@everywhere using LinearAlgebra, BallArithmetic, JLD
@everywhere D = JLD.load("../../BlaschkeMatrixSchur512.jld")

λ = D["S"].values[1]
R = 0.01

io = open("./logs/log_$(location)_Blashke_$(λ)_$(R)_$datetime.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)
@info "Added $nprocs processes"
Sys.cpu_summary(io)

@info "Schur decomposition errors"
errF = D["errF"]
errT = D["errT"]
norm_Z = D["norm_Z"]
norm_Z_inv = D["norm_Z_inv"]
@info "E_M", errF
@info "E_T", errT
@info "norm_Z", norm_Z
@info "norm_Z_inv", norm_Z_inv
N = size(D["P"])[1]

@everywhere const T_global = BallMatrix(D["S"].T)

const job_channel = RemoteChannel(()->Channel{Tuple{Int, ComplexF64}}(1024))
const result_channel = RemoteChannel(()->Channel{NamedTuple}(1024))

const certification_log = DataFrame(
    i = Int[],
    val = Ball{Float64, Float64}[],
    lo_val = Float64[],
    res = Ball{Float64, Float64}[],
    hi_res = Float64[],
    second_val = Ball{Float64, Float64}[],
    z = ComplexF64[],
    t = Float64[],
    id = Int[]
)

include("../script_functions_2.jl")

foreach(
        pid -> remote_do(dowork, pid, job_channel, result_channel),
        workers()
    )

N = 128
θs = range(0, 2π, length = N + 1)[1:(end - 1)]
zs = λ .+ R .* exp.(1im .* θs)
arcs = [(zs[i], zs[mod1(i + 1, N)]) for i in 1:N]
cache = Dict{ComplexF64, Any}()

@info "Certifying a circle of radius $R around $λ, initial partition in $N arcs"

#@info arcs
adaptive_arcs!(arcs, cache, 0.1)

function lo(x::Ball)
    lo = setrounding(Float64, RoundUp) do
            return x.c - x.r
    end
    return lo
end

JLD2.@save "./logs/certification_log_$(location)_Blashke_$(λ)_$(R)_$datetime.jld2" certification_log
CSV.write("./logs/certification_log_$(location)_Blashke_$(λ)_$(R)_$datetime.csv", certification_log)

@info "The smallest singular value along the arc is bounded below by $(minimum(certification_log.lo_val))"
l2pseudo = maximum(certification_log.hi_res)
@info "The resolvent norm for the Schur matrix in l2 norm is bounded above by $(l2pseudo)"

bound_res_original = setrounding(Float64, RoundUp) do
    
    norm_Z_sup = (norm_Z - 1).c + (norm_Z - 1).r
    norm_Z_inv_sup = (norm_Z_inv - 1).c + (norm_Z_inv - 1).r

    ϵ = max(max(errF, errT), max(norm_Z_sup, norm_Z_inv_sup))
    @info "The ϵ in the Schur theorems $ϵ"
    return (2 * (1 + ϵ^2) * l2pseudo * sqrt(N)) / (1 - 2 * ϵ * (1 + ϵ^2) * l2pseudo)
end

@info "The l1 resolvent norm for the original discretized operator is bounded above by $bound_res_original"

rmprocs(procs)