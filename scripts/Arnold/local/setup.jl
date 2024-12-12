using Logging, Dates

datetime = Dates.now()
io = open("log_local_Arnold_$datetime.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)

import Pkg
Pkg.activate("./")

using Distributed, FileIO
using FileIO, CSV
nprocs = 2
procs = addprocs(nprocs, enable_threaded_blas = true)

@info "Added $nprocs processes"
Sys.cpu_summary(io)

@everywhere import Pkg;
@everywhere Pkg.activate("./")
@everywhere Pkg.instantiate();
@everywhere Pkg.precompile()

@everywhere using JLD
@everywhere using LinearAlgebra
@everywhere ENV["OPENBLAS_NUM_THREADS"] = 4

@everywhere using BallArithmetic

@everywhere D = load("../../../ArnoldMatrixSchur256.jld")

include("../../script_functions.jl")