function compute_enclosure_arc(D, λ, ρ, r_pearl; csvfile = "", start_angle, stop_angle)
    jobs = RemoteChannel(() -> Channel{Tuple}(32))
    results = RemoteChannel(() -> Channel{NamedTuple}(32))

    Ntot = ExperimentsPseudospectra.compute_steps(
        ρ, r_pearl; arc = stop_angle - start_angle)
    @info "$Ntot svd need to be computed to certify the arc centered at $(λ), with radius $(ρ), with pearls of size $(r_pearl)"
    @info "with start angle $(start_angle) and stop angle $(stop_angle)"

    @async ExperimentsPseudospectra.submit_job(
        λ, ρ, r_pearl, jobs; start_angle = start_angle, stop_angle = stop_angle)

    foreach(
        pid -> remote_do(ExperimentsPseudospectra.dowork, pid, D["P"], jobs, results),
        workers()
    )

    avg_time = 0.0
    N = Ntot
    count = 0
    min_svd = 100
    l2pseudo = 0.0

    if csvfile == ""
        csvfile = "results_$(now())_$(λ)_$(ρ)_$(r_pearl)_$(start_angle)_$(stop_angle).csv"
    end

    tot_time = @elapsed while N > 0 # print out results
        x = take!(results)

        sv = setrounding(Float64, RoundDown) do
            return x.val_c - x.val_r
        end

        l2pseudoloc = setrounding(Float64, RoundUp) do
            return 1.0 / sv
        end

        l2pseudo = max(l2pseudo, l2pseudoloc)
        min_svd = min(min_svd, sv)

        avg_time += x.t
        N = N - 1
        count += 1

        CSV.write(csvfile, [x]; append = isfile(csvfile))

        # if count % 10000 == 0
        #     @info min_svd
        #     @info x.t, x.c
        #     @info count, ((avg_time / count) * Ntot) / (3600 * length(workers()))
        #     #@info "Done ", (1-Float64(N)/Ntot)*100, "%"
        # end
    end

    @info "The minimum singular value along the arc centered at $(λ), with radius $(ρ), with pearls of size $(r_pearl) with start angle $(start_angle) and stop angle $(stop_angle) is $(min_svd), the maximum of the l2 pseudospectra is bounded by $(l2pseudo)"
    avg_time /= Ntot
    @info "Average time for certifying an SVD", avg_time
    @info "Total time for certifying the arc", tot_time
    return min_svd, l2pseudo
end