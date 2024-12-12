function build_matrix_blaschke(K)
    r = IntervalArithmetic.interval(9)/10
    ϕ = IntervalArithmetic.interval(π)/4

    D(x) = 0.5+atan((RigorousInvariantMeasures.sinpi(2*x)-r*sin(ϕ))/(RigorousInvariantMeasures.cospi(2*x)-r*cos(ϕ)))/pi


    function T(x; c)
        return 2 * x + c * RigorousInvariantMeasures.sinpi(2 * x) +
        sqrt(IntervalArithmetic.IntervalArithmetic.interval(2)) / 2
    end
    c = 1 / (2 * IntervalArithmetic.interval(pi)) - 1 / 16
    FourierBasis = RigorousInvariantMeasures.FourierAdjoint(K, 32768)
    P = RigorousInvariantMeasures.DiscretizedOperator(FourierBasis, x -> T(x; c = c))
    return P
end

function prepare_and_output_Blaschke(K, filename = "BlaschkeMatrixSchur$K.jld")
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