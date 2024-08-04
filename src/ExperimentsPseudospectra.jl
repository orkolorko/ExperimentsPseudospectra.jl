module ExperimentsPseudospectra

import Pkg
using SetRounding

function setup()
    Pkg.instantiate()
end

# import RigorousInvariantMeasures, IntervalArithmetic

# function build_matrix_arnold(K)
#     function T(x; c)
#         return 2 * x + c * RigorousInvariantMeasures.sinpi(2 * x) +
#         sqrt(IntervalArithmetic.interval(2)) / 2
#     end
#     c = 1 / (2 * IntervalArithmetic.interval(pi)) - 1 / 16
#     FourierBasis = RigorousInvariantMeasures.FourierAdjoint(K, 32768)
#     P = RigorousInvariantMeasures.DiscretizedOperator(FourierBasis, x -> T(x; c = c))
#     return P
# end

import BallArithmetic

# function convert_matrix(P)
#     midI = IntervalArithmetic.mid
#     radI = IntervalArithmetic.radius
#     midP = midI.(real.(P.L)) + im * midI.(imag.(P.L))
#     radP = setrounding(Float64, RoundUp) do
#         return sqrt.(radI.(real.(P.L))^2 + radI.(imag.(P.L))^2)
#     end
#     return BallArithmetic.BallMatrix(midP, radP)
# end

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

export prepare_and_output
function prepare_and_output(K, filename = "ArnoldMatrixSchur$K.jld")
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

@inline compute_steps(ρ, r_pearl) = ceil(Int64, (2*pi*ρ)/r_pearl)

function submit_job(λ, ρ, r_pearl, job_queue; N = compute_steps(ρ, r_pearl), start = 0, stop = N)
    @info "$(stop-start) jobs"
    for i in start:stop
        put!(job_queue, (i, λ+ρ*exp(2*pi*im*i/N), r_pearl))
    end
    return N
end

import Distributed

function dowork(P, jobs, results)
    while true
        i, c, r_pearl = take!(jobs)
        @info c, r_pearl
        z = BallArithmetic.Ball(c, r_pearl)
        t = @elapsed Σ = BallArithmetic.svdbox(P-z*LinearAlgebra.I)
        put!(results, (i = i, val = Σ[end], second_val = Σ[end-1] , c = c, r_pearl = r_pearl, t = t, id = Distributed.myid()))
    end
end

end
