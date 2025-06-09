bound_res_original = setrounding(Float64, RoundUp) do
    
    norm_Z_sup = (norm_Z - 1).c + (norm_Z - 1).r
    norm_Z_inv_sup = (norm_Z_inv - 1).c + (norm_Z_inv - 1).r

    ϵ = max(max(errF, errT), max(norm_Z_sup, norm_Z_inv_sup))
    return (2 * (1 + ϵ^2) * l2pseudo * sqrt(N)) / (1 - 2 * ϵ * (1 + ϵ^2) * l2pseudo)
end

@info "The norm of the resolvent for the discretized operator is", bound_res_original

rmprocs(procs)
close(io)