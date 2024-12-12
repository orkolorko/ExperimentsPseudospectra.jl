module ExperimentsPseudospectra

import Pkg
using SetRounding

function setup()
    Pkg.instantiate()
end

import RigorousInvariantMeasures, IntervalArithmetic
import BallArithmetic

using LinearAlgebra
using JLD


 
include("utilities.jl")

include("parallel.jl")

include("Arnold.jl")
include("Blaschke.jl")
export prepare_and_output_Arnold, prepare_and_output_Blaschke



end
