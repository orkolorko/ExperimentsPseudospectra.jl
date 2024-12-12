@inline compute_steps(ρ, r_pearl; arc = 2*pi) = ceil(Int64, (arc * ρ) / r_pearl)

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