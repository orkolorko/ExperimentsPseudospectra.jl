function build_matrix_arnold(K)
    function T(x; c)
        return 2 * x + c * RigorousInvariantMeasures.sinpi(2 * x) + 0.25
    end
    b = 5/64+1/128+1/256
    c = (1 / (2 * IntervalArithmetic.interval(pi)) - b)
    FourierBasis = RigorousInvariantMeasures.FourierAdjoint(K, 1048576)
    P = RigorousInvariantMeasures.DiscretizedOperator(FourierBasis, x -> T(x; c = c))
    return P
end

function prepare_and_output_Arnold(K, filename = "ArnoldMatrixSchur$K.jld2")
    P = ExperimentsPseudospectra.convert_matrix(ExperimentsPseudospectra.build_matrix_arnold(K))
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