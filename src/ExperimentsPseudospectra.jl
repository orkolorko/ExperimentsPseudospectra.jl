module ExperimentsPseudospectra

import Pkg
using SetRounding

function setup()
    Pkg.instantiate()
end

import RigorousInvariantMeasures, IntervalArithmetic

function build_matrix_arnold(K)
     function T(x; c)
         return 2 * x + c * RigorousInvariantMeasures.sinpi(2 * x) +
         sqrt(IntervalArithmetic.interval(2)) / 2
     end
     c = 1 / (2 * IntervalArithmetic.interval(pi)) - 1 / 16
     FourierBasis = RigorousInvariantMeasures.FourierAdjoint(K, 32768)
     P = RigorousInvariantMeasures.DiscretizedOperator(FourierBasis, x -> T(x; c = c))
     return P
end

function build_matrix_Blaschke(K)
    r = interval(9)/10
    ϕ = interval(π)/4

    D(x) = 0.5+atan((RigorousInvariantMeasures.sinpi(2*x)-r*sin(ϕ))/(RigorousInvariantMeasures.cospi(2*x)-r*cos(ϕ)))/pi


    function T(x; c)
        return 2 * x + c * RigorousInvariantMeasures.sinpi(2 * x) +
        sqrt(IntervalArithmetic.interval(2)) / 2
    end
    c = 1 / (2 * IntervalArithmetic.interval(pi)) - 1 / 16
    FourierBasis = RigorousInvariantMeasures.FourierAdjoint(K, 32768)
    P = RigorousInvariantMeasures.DiscretizedOperator(FourierBasis, x -> T(x; c = c))
    return P
end

import BallArithmetic

function convert_matrix(P)
    midI = IntervalArithmetic.mid
    radI = IntervalArithmetic.radius
    midP = midI.(real.(P.L)) + im * midI.(imag.(P.L))
    radP = setrounding(Float64, RoundUp) do
        return sqrt.(radI.(real.(P.L))^2 + radI.(imag.(P.L))^2)
    end
    return BallArithmetic.BallMatrix(midP, radP)
end

using LinearAlgebra

function compute_schur_and_error(A::BallArithmetic.BallMatrix)
    S = LinearAlgebra.schur(Complex{Float64}.(A.c))

    bZ = BallArithmetic.BallMatrix(S.Z)
    errF = BallArithmetic.svd_bound_L2_opnorm(bZ' * bZ - I)

    bT = BallArithmetic.BallMatrix(S.T)
    errT = BallArithmetic.svd_bound_L2_opnorm(bZ * bT * bZ' - A)

    sigma_Z = BallArithmetic.svdbox(bZ)

    norm_Z = sigma_Z[1]
    norm_Z_inv = 1.0 / sigma_Z[end]

    return S, errF, errT, norm_Z, norm_Z_inv
end

using JLD

export prepare_and_output_Arnold, prepare_and_output_Blaschke
 
function prepare_and_output_Arnold(K, filename = "ArnoldMatrixSchur$K.jld")
    P = ExperimentsPseudospectra.convert_matrix(ExperimentsPseudospectra.build_matrix_arnold(256))
    S, errF, errT, norm_Z, norm_Z_inv = compute_schur_and_error(P)
    jldopen(filename, "w") do file
        write(file, "P", P)
        write(file, "S", S)
        write(file, "errF", errF)
        write(file, "errT", errT)
        write(file, "norm_Z", norm_Z)
        write(file, "norm_Z_inv", norm_Z_inv)
    end
end

@inline compute_steps(ρ, r_pearl; arc = 2*pi) = ceil(Int64, (arc * ρ) / r_pearl)

# function submit_job(
#         λ, ρ, r_pearl, job_queue; N = compute_steps(ρ, r_pearl), start = 0, stop = N)
#     @info "$(stop-start) jobs"
#     for i in start:stop
#         put!(job_queue, (i, λ + ρ * exp(2 * pi * im * i / N), r_pearl))
#     end
#     return N
# end

function submit_job(
    λ, ρ, r_pearl, job_queue; start_angle = 0, stop_angle = 2*pi, N = compute_steps(ρ, r_pearl; arc = stop_angle-start_angle))
    for (i, θ) in enumerate(range(; start = start_angle, stop = stop_angle, length = N))
        # cis(x)
        # More efficient method for exp(im*x) by using Euler's formula: cos(x) + i sin(x) = \exp(i x)
        put!(job_queue, (i, θ,  λ + ρ * cis(θ), r_pearl))
    end
    return N
end


import Distributed

function dowork(P, jobs, results)
    while true
        i,θ, c, r_pearl = take!(jobs)
        #@info c, r_pearl
        z = BallArithmetic.Ball(c, r_pearl)
        t = @elapsed Σ = BallArithmetic.svdbox(P - z * LinearAlgebra.I)
        put!(results,
            (i = i,
                val_c = Σ[end].c,
                val_r = Σ[end].r,
                second_val_c = Σ[end - 1].c,
                second_val_r = Σ[end - 1].r,
                c = c,
                radian = θ,
                r_pearl = r_pearl,
                t = t,
                id = Distributed.myid()))
    end
end

end
