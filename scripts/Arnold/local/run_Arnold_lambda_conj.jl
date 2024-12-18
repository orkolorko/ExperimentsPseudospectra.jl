include("setup.jl")

filename = "Arnold_lambda_conj"
filename*= "$datetime"


S = D["S"]
ρ = 0.001
r_pearl = (4.9 * 10^(-9))
λ = S.values[2]

@info "Certifying ", λ, "radius", ρ, "radius pearl", r_pearl

start_angle = 0
stop_angle = 2*pi

compute_enclosure_arc(
        D, λ, ρ, r_pearl; start_angle = start_angle, stop_angle = stop_angle, csvfile = "$filename.csv")

save("$filename.jld", "lambda", λ, "rho", ρ, "r_pearl", r_pearl, "min_svd", min_svd, "l2pseudo", l2pseudo)

include("cleaning.jl")
