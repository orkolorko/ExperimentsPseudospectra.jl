include("setup.jl")

filename = filename = "Arnold_one"
filename*= "$datetime"

S = D["S"]
ρ = 0.1
r_pearl = (4.9 * 10^(-2))
λ = 1.0

@info "Certifying ", λ, "radius", ρ, "radius pearl", r_pearl

start_angle = 0
stop_angle = 2 * pi

min_svd, l2pseudo = compute_enclosure_arc(
    D, λ, ρ, r_pearl; start_angle = start_angle,
    stop_angle = stop_angle, csvfile = "$filename.csv")

save("$filename.jld", "lambda", λ, "rho", ρ, "r_pearl", r_pearl, "min_svd", min_svd, "l2pseudo", l2pseudo)

include("cleaning.jl")