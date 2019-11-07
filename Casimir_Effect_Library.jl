using QuadGK
using Plots
using LaTeXStrings
using LinearAlgebra

## Parameters
const ν = 1e-3;         # Relative tolerance for integration
const η = 1e-7;         # Small number used for i0
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

# Impurity type
struct Impurity
    pos :: Int        # Position of the impurity
    M   :: Float64    # Impurity mass in units of the chain atom mass
end

# Function for producing an array of impurities
function impurity_cluster(n :: Int, pos, M)
    imps = map(x -> Impurity(pos + x - 1, M), 1 : n)
    return imps
end

## Functions
# Bose-Einstein distribution
function nB(x :: Float64)
    return (1 ./ (exp(x) - 1))
end

# Propagator for the finite-T energy
function P(D :: Int, z :: Float64)
    x = (z + 1im * η)^2
    res = -(D == 0) + √(-x / (1 - x)) * (1 - 2 * x - 2 * √(-x) * √(1 - x))^D
    return res
end

# Propagator for the zero-T energy
function P_T0(D :: Int, θ :: Float64)
    res = -(D == 0) +  tanh(θ) * exp(-2 * D * θ)
    return res
end

# Scattering matrix for finite-T energy
function Δ(z, Imps :: Array{Impurity, 1})
    nImps = length(Imps)
    Ms = map(x -> x.M, Imps)
    α = Diagonal(1 .- 1 ./ Ms)

    Imp_Mat = repeat(Imps, 1, nImps)    # Impurity position matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity position matrix

    Prop_Mat = map((x, y) -> P(abs.(x.pos - y.pos), z), Imp_Mat, ImpT_Mat)

    return (Matrix{Int}(I, nImps, nImps) .+ Prop_Mat * α)
end


# Scattering matrix for zero-T energy
function Δ_T0(θ, Imps :: Array{Impurity, 1})
    nImps = length(Imps)
    Ms = map(x -> x.M, Imps)
    α = Diagonal(1 .- 1 ./ Ms)

    Imp_Mat = repeat(Imps, 1, nImps)    # Impurity position matrix
    ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity position matrix

    Prop_Mat = map((x, y) -> P_T0(abs.(x.pos - y.pos), θ), Imp_Mat, ImpT_Mat)

    return (Matrix{Int}(I, nImps, nImps) .+ Prop_Mat * α)
end


# FI Integrand for finite-T energy
function F_I_Integrand(Imps, T, z)
    Ms = map(x -> x.M, Imps)
    α = 1 .- 1 ./ Ms

    Δ0_Inv = Diagonal(1 ./ (1 .+ P(0, z) .* α))
    Δ_ = Δ(z, Imps)

    return (imag(log(Complex(det(Δ0_Inv * Δ_)))) * (nB(z / T) ))
end

# FI Integrand for zero-T energy
function F_I_Integrand_T0(Imps, θ)
    Ms = map(x -> x.M, Imps)
    α = 1 .- 1 ./ Ms

    Δ0_Inv = Diagonal(1 ./ (1 .+ P_T0(0, θ) .* α))
    Δ_ = Δ_T0(θ, Imps)

    return real(log(Complex(det(Δ0_Inv * Δ_))) * cosh(θ))
end

# FI for finite T
function F_I(Imps, T)
    f(z) = F_I_Integrand(Imps, T, z)
    res = quadgk(f, -Inf, 0, Inf, maxevals = NumEvals)[1]
    return (res / (2 * π))
end

# FI for zero T
function F_I_T0(Imps)
    f(θ) = F_I_Integrand_T0(Imps, θ)
    res = quadgk(f, 0, 1e2, maxevals = NumEvals)[1]
    return (res / (2 * π))
end

# Interaction energy for a pair of clusters of length n and mass M.
# The separation D is the distance between the rightmost mass of the left
# cluster and the leftmost mass of the right one.
# THIS IS COMPUTED AT T = 0
function F_I_cluster(D, n, M)
    left_cluster = impurity_cluster(n, 0, M)
    right_cluster = impurity_cluster(n, n - 1 + D, M)
    imps = vcat(left_cluster, right_cluster)
    res = F_I_T0(imps)
    return res
end

# Free Energy from Exact Diagonalization for two impurities
function Exact_Free_Energy(N, M1, M2, D, T)
    # Prepare a pristine chain potential energy matrix
    dv = 2 .* ones(N)
    ev = -ones(N - 1)

    U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1,  N] = -1
    U_Mat[N,  1] = -1

    # Prepare the matrix of masses
    M_Mat = Diagonal(ones(N))
    # Replace two pristine masses by impurities
    M_Mat[1 + D, 1 + D] = M1
    M_Mat[1, 1] = M2

    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = real(eigvals(inv(M_Mat) * U_Mat))
    # These Ω's are given in units of sqrt(K / m). We, however, want to give
    # the energies in the units of 2*sqrt(K / m). Hence, we divide the
    # energies by 2:
    Ω = real(sqrt.(Ω2[2 : N]) / 2)

    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))

    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2

    return total_energy
end

# Colors for plotting
my_red = RGB(215/255,67/255,84/255)
my_green = RGB(106/255,178/255,71/255)
my_blue = RGB(100/255,101/255,218/255)
my_violet = RGB(169/255,89/255,201/255)
my_orange = RGB(209/255,135/255,46/255)

colors = [my_red
        , my_green
        , my_blue
        , my_violet
        , my_orange
         ]
