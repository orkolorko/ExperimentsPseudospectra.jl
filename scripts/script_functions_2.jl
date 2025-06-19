@everywhere function dowork(jobs, results)
    while true
        i, z = take!(jobs)
        @debug "Received and working on," z
        bz = BallArithmetic.Ball(z, 0.0)

        t = @elapsed Σ = BallArithmetic.svdbox(T_global - bz * LinearAlgebra.I)

        val = Σ[end]
        res = 1 / val

        lo_val = setrounding(Float64, RoundDown) do
            return val.c - val.r
        end

        hi_res = setrounding(Float64, RoundUp) do
            return res.c + res.r
        end

        put!(results,
            (
                i = i,
                val = val,
                lo_val = lo_val,
                res = res,
                hi_res = hi_res,
                second_val = Σ[end - 1],
                z = z,
                t = t,
                id = myid()
            ))
    end
end

# --- Adaptive method ---
function adaptive_arcs!(arcs::Vector{Tuple{ComplexF64, ComplexF64}},
        cache::Dict{ComplexF64, Any},
        pending::Dict{Int, Tuple{ComplexF64, ComplexF64}},
        id_counter,
        η::Float64;
        check_interval = 1000)
    i = 1
    processed = 0
    new = 0

    cycle = true
    @info "Running, initial length arcs", length(arcs)

    while !isempty(arcs)
        z_a, z_b = arcs[i]
        if haskey(cache, z_a)
            σ_a = cache[z_a][1]
        else
            job_id = id_counter
            @debug "put on job channel", z_a, z_b
            put!(job_channel, (job_id, z_a))
            pending[job_id] = (z_a, z_b)
            id_counter += 1
            deleteat!(arcs, i)
            continue
        end

        ℓ = abs(z_b - z_a)
        ε = ℓ / σ_a

        sup_ε = setrounding(Float64, RoundUp) do
            return ε.c + ε.r
        end

        if sup_ε <= η
            deleteat!(arcs, i)
        else
            z_m = (z_a + z_b) / 2
            deleteat!(arcs, i)
            insert!(arcs, i, (z_m, z_b))
            insert!(arcs, i, (z_a, z_m))
            new += 1
        end

        processed += 1
        if processed % check_interval == 0
            @info "Processed $processed arcs..."
            @info "Length of arcs to be processed", length(arcs)
            @info "Pending", length(pending)
            @info "Proportion of new arcs", new / check_interval
            new = 0

            flush(io)
            @debug "Processed", processed
            @debug "Waiting for pending"
            while isready(result_channel)
                result = take!(result_channel)
                z = result.z
                cache[z] = (result.val, result.second_val)
                z_a, z_b = pending[result.i]
                delete!(pending, result.i)
                push!(arcs, (z_a, z_b))
                push!(certification_log, result)
            end
            cycle = save_snapshot!(arcs, cache, certification_log, pending, snapshot, cycle)
        end
    end

    # Final drain of remaining results
    while !isempty(pending)
        result = take!(result_channel)
        z = result.z
        cache[z] = (result.val, result.second_val)
        z_a, z_b = pending[result.i]
        delete!(pending, result.i)
        push!(arcs, (z_a, z_b))
        push!(certification_log, result)
    end

    @debug "After initial process", length(arcs)
    # Recursively refine new arcs
    if !isempty(arcs)
        adaptive_arcs!(arcs, cache, pending, id_counter, η; check_interval)
    end
end

function bound_res_original(l2pseudo, η, norm_Z, norm_Z_inv, errF, errT, N)
    bound = setrounding(Float64, RoundUp) do
        l2pseudo = l2pseudo * 1 / (1 - η)
        norm_Z_sup = (norm_Z - 1).c + (norm_Z - 1).r
        norm_Z_inv_sup = (norm_Z_inv - 1).c + (norm_Z_inv - 1).r

        ϵ = max(max(errF, errT), max(norm_Z_sup, norm_Z_inv_sup))
        @info "The ϵ in the Schur theorems $ϵ"
        return (2 * (1 + ϵ^2) * l2pseudo * sqrt(N)) / (1 - 2 * ϵ * (1 + ϵ^2) * l2pseudo)
    end
    return bound
end

function choose_snapshot_to_load(basepath::String)
    files = [basepath * "_A.jld2", basepath * "_B.jld2"]
    valid_files = filter(isfile, files)
    if isempty(valid_files)
        return nothing  # No prior snapshot found
    end
    # Sort by modification time descending
    sorted = sort(valid_files, by = f -> stat(f).mtime, rev = true)
    try
        snapshot = JLD2.load(sorted[1])
        return snapshot
    catch e
        @warn "Could not load snapshot file $(sorted[1]), possibly corrupted. Trying backup." exception=(
            e, catch_backtrace())
        if length(sorted) > 1
            try
                snapshot = JLD2.load(sorted[2])
                return snapshot
            catch e
                @warn "Both snapshot files failed to load." exception=(e, catch_backtrace())
                return nothing
            end
        else
            return nothing
        end
    end
end

function save_snapshot!(arcs, cache, log, pending, basepath::String, toggle::Bool)
    filename = basepath * (toggle ? "_A.jld2" : "_B.jld2")
    JLD2.@save filename arcs cache log pending
    return !toggle  # Flip toggle
end