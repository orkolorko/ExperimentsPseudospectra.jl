using Logging, Dates

datetime = Dates.now()
io = open("log_local_Arnold_$datetime.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)

import Pkg
Pkg.activate("./")

using Distributed, FileIO
using FileIO, CSV
nprocs = 8

procs = addprocs(nprocs, enable_threaded_blas = false)

@info "Added $nprocs processes"
Sys.cpu_summary(io)

@everywhere import Pkg;
@everywhere Pkg.activate("./")
@everywhere Pkg.instantiate();
@everywhere Pkg.precompile()

@everywhere using JLD
@everywhere using LinearAlgebra
@everywhere ENV["OPENBLAS_NUM_THREADS"] = 1

@everywhere using BallArithmetic

@everywhere D = load("../../../ArnoldMatrixSchur256.jld")

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

include("../../script_functions.jl")