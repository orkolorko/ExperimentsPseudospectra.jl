include("setup.jl")

filename = "Arnold_center"

S = D["S"]
λ = 0.0
ρ = 0.27812055400531616

eig1 = S.values[1]
eig2 = S.values[2]

# we identify two sectors near the eigenvalues where we are going to work
# with a smaller step 

size_sec = 0.0001
angle_1 = angle(eig1)
angle_2 = 2 * pi + angle(eig2)

sectors = [(0.0, angle_1 - size_sec, 2^(-23))
           (angle_1 - size_sec, angle_1 + size_sec, 2^(-24));
           (angle_1 + size_sec, angle_2 - size_sec, 2^(-23));
           (angle_2 - size_sec, angle_2 + size_sec, 2^(-24));
           (angle_2 + size_sec, 2*pi, 2^(-23))]


out = []

for S in sectors
    start = S[1]
    stop = S[2]
    r_pearl = S[3]
    min_svd, l2pseudo = compute_enclosure_arc(
        D, λ, ρ, r_pearl; start_angle = π / 2, stop_angle = angle_1 - size_sec, csvfile = "$filename.csv")
    push!(out, (min_svd, l2pseudo))
end


save("$filename.jld", "lambda", λ, "rho", ρ, "sectors", sectors, "output", out)

include("cleaning.jl")