using IntervalArithmetic

function build_matrix_blaschke(K)
    r = sqrt(interval(2))*(interval(3)/8) 
    ϕ = interval(π) / 8

    D(x) = 0.5+atan((RigorousInvariantMeasures.sinpi(2*x)-r*sin(ϕ))/(RigorousInvariantMeasures.cospi(2*x)-r*cos(ϕ)))/pi

    FourierBasis = RigorousInvariantMeasures.FourierAdjoint(K, 1048576)
    P = RigorousInvariantMeasures.DiscretizedOperator(FourierBasis, x -> D(x))
    return P
end

function prepare_and_output_Blaschke(K, filename = "BlaschkeMatrixSchur$K.jld2")
    P = ExperimentsPseudospectra.convert_matrix(ExperimentsPseudospectra.build_matrix_blaschke(K))
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