# using LaTeXStrings
# using LinearAlgebra
# using Plots
using QuadGK

## Parameters
const ν = 1e-4;         # Small number for relative tolerance
const α = 1e-8;         # Small number for absolute tolerance
const η = 1e-5;         # Small number for moving the contour off the real axis
const NumEvals = 1e7;   # Maximum number of evaluations in quadgk

# Impurity type
struct Impurity
    pos::Int        # Position of the impurity
    M::Float64      # Impurity mass in units of m
    Δ::Float64      # Harmonic potential perturbation in units of k
end

function Ω(θ, K)
    return √(4 * K * sin(θ / 2)^2 + 1)
end


function Π_jl(z, K, D)

    f_int(θ) = 1 ./ (-1 .+ (Ω(θ, K) ./ z)^2) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π)[1]
    return (res / (2 * π * z^2))
end


@time Π_jl(.09im, 9, 8)

#
#
# function E_I_Integrand(ω, system)
#     Ms = system.Ms
#     imps = system.Imps
#     K = system.K
#     imp_host_idx = map(x -> x.n, imps)
#     m = map(x -> Ms[x], imp_host_idx) |> Diagonal
#     m_sqrt = sqrt(m)
#     m_sqrt_inv = inv(m_sqrt)
#
#     nImps = length(imps)
#     imps_mat = repeat(imps, 1, nImps)
#     imps_mat_T = permutedims(imps_mat)
#
#     Π = map((x, y) -> Π_jl(ω, Ms, K, x, y), imps_mat, imps_mat_T)
#     # Matrix to remove couplings between defects
#     mask = Diagonal(ones(nImps))
#     ΠD = mask .* Π
#
#     Ξ = [
#         m_sqrt_inv*Π*m_sqrt_inv ω*m_sqrt_inv*Π*m_sqrt
#         ω*m_sqrt*Π*m_sqrt_inv (ω^2)*m_sqrt*Π*m_sqrt+m
#     ]
#     ΞD = [
#         m_sqrt_inv*ΠD*m_sqrt_inv ω*m_sqrt_inv*ΠD*m_sqrt
#         ω*m_sqrt*ΠD*m_sqrt_inv (ω^2)*m_sqrt*ΠD*m_sqrt+m
#     ]
#
#     unit_mat = ones(2 * length(imps)) |> Diagonal
#     Δ_Λ =
#         vcat(
#             map(x -> x.Δ, imps),
#             map(x -> 1 / x.M - 1 / system.Ms[x.n], imps),
#         ) |> Diagonal
#     return (
#         (det(unit_mat .+ Ξ * Δ_Λ) |> Complex |> log) -
#         (det(unit_mat .+ ΞD * Δ_Λ) |> Complex |> log)
#     )
# end



#
# function propagator(z, D)
#     x = z^2
#     res = ((1 - 2 * x - 2 * √(-x) * √(1 - x))^D) / (√(-x) * √(1 - x))
#     return res
# end
#
# function F_I_Integrand(z, imps)
#     nImps = length(imps)
#     Λ_ = 1 ./ map(x -> x.M, imps) .- 1 |> Diagonal
#     Δ_ = map(x -> x.Δ, imps) |> Diagonal
#
#     Imp_Mat = repeat(imps, 1, nImps)    # Impurity matrix
#     ImpT_Mat = permutedims(Imp_Mat)     # Transpose of impurity matrix
#
#     Y = map((x, y) -> propagator(z, abs.(x.pos - y.pos)), Imp_Mat, ImpT_Mat)
#
#     return (
#                Matrix{Int}(I, 2 * nImps, 2 * nImps) + [
#                    (Y * Δ_) (z .* Y * Λ_)
#                    (z .* Y * Δ_) (
#                        (Matrix{Int}(I, nImps, nImps) + z^2 .* Y) * Λ_
#                    )
#                ]
#            ) |>
#            det |>
#            log
# end
#
# function F_I(imps)
#     f(x) = F_I_Integrand(1im * x, imps)
#     res = quadgk(f, 0, Inf, maxevals = NumEvals)[1]
#     return (res / (2 * π))
# end
#
# function Exact_Free_Energy(N, imps)
#     # Prepare a pristine chain potential energy matrix
#     dv = 2 .* ones(N)
#     ev = -ones(N - 1)
#
#     U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
#     # Make the system periodic by coupling the first and the last mass
#     U_Mat[1, N] = -1
#     U_Mat[N, 1] = -1
#     # Add the external harmonic potential. The factor of 4 is because Δ is
#     # given in units of 4 * K
#     for ii = 1:length(imps)
#         imp = imps[ii]
#         U_Mat[imp.pos+1, imp.pos+1] += 4 * imp.Δ
#     end
#
#     # Prepare the matrix of masses
#     M_Mat = Diagonal(ones(N))
#     # Modify the masses
#     for ii = 1:length(imps)
#         imp = imps[ii]
#         M_Mat[imp.pos+1, imp.pos+1] = imp.M
#     end
#
#     # Calculate the eigenvalues of the system which give
#     # the squares of the energies
#     Ω2 = real(eigvals(inv(M_Mat) * U_Mat))
#     # These Ω's are given in units of sqrt(K / m). We, however, want to give
#     # the energies in the units of 2*sqrt(K / m). Hence, we divide the
#     # energies by 2:
#     Ω = real(sqrt.(Ω2[2:N]) / 2)
#
#     return (sum(Ω) / 2)
# end
#
#
# ## Colors for plotting
# my_yellow = RGB(255 / 255, 255 / 255, 153 / 255)
# my_green = RGB(127 / 255, 201 / 255, 127 / 255)
# my_blue = RGB(56 / 255, 108 / 255, 176 / 255)
# my_violet = RGB(190 / 255, 174 / 255, 212 / 255)
# my_red = RGB(240 / 255, 2 / 255, 127 / 255)
# my_orange = RGB(253 / 255, 192 / 255, 134 / 255)
#
# colors = [my_red, my_green, my_blue, my_violet, my_orange]
