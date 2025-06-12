@everywhere function dowork(jobs, results)
    while true
        i, z = take!(jobs)
        @debug "Received and working on," z
        bz = BallArithmetic.Ball(z, 0.0)

        t = @elapsed Σ = BallArithmetic.svdbox(T_global - bz * LinearAlgebra.I)

        val = Σ[end]
        res = 1/val

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
        η::Float64)

    pending = Dict{Int, Tuple{ComplexF64, ComplexF64}}()
    id_counter = 1
    i = 1
    processed = 0
    check_interval = 100

    toggle = true
    @debug "Running, initial length arcs", length(arcs)

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
        end

        processed += 1
        if processed % check_interval == 0
            @info "Processed $processed arcs..."
            @info "Length of arcs to be processed", length(arcs)
            @info "Pending", length(pending)
            
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
            if toggle
                filename_processed = filename*"_A.jld2"
                toggle = false
            else
                filename_processed = filename*"_B.jld2"
                toggle = true
            end

            JLD2.@save filename_processed certification_log
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
        adaptive_arcs!(arcs, cache, η)
    end
end