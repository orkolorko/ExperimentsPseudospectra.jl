#!/usr/bin/env julia
# Regenerate Figures 1, 3, 4 with equal aspect ratio (axis equal).
# Uses precomputed Schur decompositions from JLD2 files to avoid recomputing.
# Saves PDFs directly to the paper folder.

import Pkg
Pkg.activate(joinpath(@__DIR__, "notebook"))
Pkg.instantiate()
Pkg.add("JLD2")  # add if not present

using Plots, LinearAlgebra, LaTeXStrings, JLD2

const PAPER_DIR = "/home/isaia/Dropbox/Lavoro/Collaborators/Blumenthal-Nisoli-Taylor-Crush/DraftPseudospectra/resubmission 2-2026/resubmission with color"

circ = [cis(θ) for θ in 0:0.005:2π]   # fine circle for smooth curves

# ============================================================
# ARNOLD (perturbed doubling map)
# ============================================================
println("=== Arnold: loading precomputed Schur from ArnoldMatrixSchur128.jld2 ===")
F_arnold = JLD2.load(joinpath(@__DIR__, "ArnoldMatrixSchur128.jld2"), "S")
eigs_a = diag(F_arnold.T)

λ2_a = eigs_a[1]   # ≈ -0.052864 - 0.206250i
λ3_a = eigs_a[2]   # ≈ -0.052864 + 0.206250i
λ1_a = F_arnold.values[end]   # 1.0
println("Arnold eigenvalues: λ2=$(round(λ2_a, digits=6)), λ3=$(round(λ3_a, digits=6))")

# Figure 4b: circles (Arnold) — full view
println("Plotting circles_Arnold.pdf ...")
plot(real.(circ), imag.(circ);
     label = "R=1", linewidth=1, aspect_ratio = :equal,
     xlabel = "Re(λ)", ylabel = "Im(λ)", title = "Arnold — enclosing circles")
plot!(real.(λ1_a .+ 0.1 .* circ), imag.(λ1_a .+ 0.1 .* circ);
      label = L"F_1\ (\lambda=1)", linewidth=2)
plot!(real.(λ2_a .+ 0.001 .* circ), imag.(λ2_a .+ 0.001 .* circ);
      label = L"F_2", linewidth=2)
plot!(real.(λ3_a .+ 0.001 .* circ), imag.(λ3_a .+ 0.001 .* circ);
      label = L"F_3", linewidth=2)
plot!(real.(0.21 .* circ), imag.(0.21 .* circ);
      label = L"F_0\ (r=0.21)", linewidth=1, linestyle=:dash)
savefig(joinpath(PAPER_DIR, "circles_Arnold.pdf"))
println("  saved.")

# Figure 1: zoom around F₂/F₃ (circles_Arnold_zoom)
println("Plotting circles_Arnold_zoom.pdf ...")
plot!(xlims = (-0.1, 0.0), ylims = (0.1, 0.22), aspect_ratio = :equal)
savefig(joinpath(PAPER_DIR, "circles_Arnold_zoom.pdf"))
println("  saved.")

# ============================================================
# BLASCHKE product
# ============================================================
println("\n=== Blaschke: loading precomputed Schur from BlaschkeMatrixSchur128.jld2 ===")
F_blaschke = JLD2.load(joinpath(@__DIR__, "BlaschkeMatrixSchur128.jld2"), "S")
eigs_b = diag(F_blaschke.T)

λ2_b = eigs_b[1]   # ≈ 0.4899611 + 0.2029485i
λ3_b = eigs_b[2]   # ≈ 0.4899611 - 0.2029485i
λ1_b = F_blaschke.values[end]   # 1.0
println("Blaschke eigenvalues: λ2=$(round(λ2_b, digits=6)), λ3=$(round(λ3_b, digits=6))")

# Figure 3b: circles (Blaschke)
println("Plotting circles_Blashke.pdf ...")
plot(real.(circ), imag.(circ);
     label = "R=1", linewidth=1, aspect_ratio = :equal,
     xlabel = "Re(λ)", ylabel = "Im(λ)", title = "Blaschke — enclosing circles")
plot!(real.(λ1_b .+ 0.1 .* circ),  imag.(λ1_b .+ 0.1 .* circ);
      label = L"F_1\ (\lambda=1)", linewidth=2)
plot!(real.(λ2_b .+ 0.01 .* circ), imag.(λ2_b .+ 0.01 .* circ);
      label = L"F_2", linewidth=2)
plot!(real.(λ3_b .+ 0.01 .* circ), imag.(λ3_b .+ 0.01 .* circ);
      label = L"F_3", linewidth=2)
plot!(real.(0.51 .* circ), imag.(0.51 .* circ);
      label = L"F_0\ (r=0.51)", linewidth=1, linestyle=:dash)
savefig(joinpath(PAPER_DIR, "circles_Blashke.pdf"))
println("  saved.")

println("\nAll figures saved to:\n  $PAPER_DIR")
