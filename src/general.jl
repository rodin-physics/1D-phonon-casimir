using CairoMakie
using DelimitedFiles
using LinearAlgebra
using ProgressMeter
using QuadGK
using Statistics

## Parameters
const ν = 1e-5;             # Relative tolerance for integration
const η = 1e-5;             # Small number used for i0
const α = 1e-9;             # Small number for absolute tolerance
const max_omega = 100;    # Large number used to truncate the finite
# temperature Matsubara sum
const NumEvals = 1e7;       # Maximum number of evaluations in quadgk

# Bose-Einstein distribution
function nB(x::Float64)
    return (1 ./ (exp(x) - 1))
end

## Types
# Impurity type
struct Impurity
    pos::Int    # Position of the impurity unit cell
    M::Float64  # Impurity mass in units of m
    Δ::Float64  # External harmonic potential in units of k
end

# System type
struct System
    Imps::Vector{Impurity}  # Impurities in the system
    T::Float64              # Temperature of the system
    N::Int                  # Number of unit cells in the system (used for ED)
    K::Float64              # Confining potential for the atoms in the chain
end

# Mode energies in units of √(k / m). The confining potential K is in units of k
@inline function Ω(θ, K)
    return √(4 * sin(θ / 2)^2 + K)
end

# Π(z) for z ≠ 0
@inline function Π_jl(z, K, D)
    f_int(θ) = 1 ./ (-1 .+ (Ω(θ, K) ./ z)^2) * exp(1im * θ * D)
    res = quadgk(f_int, 0, 2 * π, atol = α)[1]
    return (res / (2 * π * z^2))
end

# Π(0)
@inline function Π0_jl(K, D)
    f_int(θ) = 1 ./ (Ω(θ, K)^2) * exp(1im * θ * D)
    res = quadgk(f_int, 0, 2 * π, atol = α)[1]
    return (res / (2 * π))
end
