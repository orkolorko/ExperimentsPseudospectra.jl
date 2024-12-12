function convert_matrix(P)
    midI = IntervalArithmetic.mid
    radI = IntervalArithmetic.radius
    midP = midI.(real.(P.L)) + im * midI.(imag.(P.L))
    radP = setrounding(Float64, RoundUp) do
        return sqrt.(radI.(real.(P.L))^2 + radI.(imag.(P.L))^2)
    end
    return BallArithmetic.BallMatrix(midP, radP)
end

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