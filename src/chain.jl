include("general.jl")

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

## Exact Diagonalization 1D to obtain the total Free Energy
function exact_F(system)
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    # Prepare the matrix of masses
    M_List = ones(N)
    # Prepare a list of unperturbed springs and a confining potential
    dv = (2 + K) .* ones(N)
    # Replace the pristine masses by the impurities
    # and add a harmonic confinement term
    for ii in Imps
        coord = ii.pos
        M_List[coord] = ii.M
        dv[coord] = dv[coord] + ii.Δ
    end

    # Assemble the matrices
    ev = -ones(N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, N] = -1
    U_Mat[N, 1] = -1
    M_Mat = Diagonal(M_List)
    # Calculate the eigenvalues of the system which give
    # the squares of the energies
    Ω2 = eigvals(inv(M_Mat) * U_Mat) |> real
    Ω = sqrt.(abs.(Ω2))  # in units √(k / μ).
    # The free energy for each mode consists of the vacuum portion Ω / 2 and
    # the finite-T portion T * log(1 - exp(-Ω / T))
    total_energy = T * sum(log.(1 .- exp.(-Ω ./ T))) + sum(Ω) / 2
    return (total_energy)
end

# -TS term for the exact diagonalization
function exact_neg_TS(system)
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K

    # Prepare a list of unperturbed springs
    dv = (2 + K) .* ones(N)
    # Prepare a list of Δs
    Δs = zeros(N)

    for ii in Imps
        Δs[ii.n] = ii.Δ
    end
    Δ_mat = Diagonal(Δs)
    # Assemble the matrices
    ev = -ones(N - 1)
    U_Mat = SymTridiagonal(dv, ev) + zeros(N, N)
    # Make the system periodic by coupling the first and the last mass
    U_Mat[1, N] = -1
    U_Mat[N, 1] = -1

    res = T * log(det(Diagonal(ones(N)) + inv(U_Mat) * Δ_mat)) / 2

    return res
end

function exact_F_I(system)
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    pristine_system = exact_F(System([], T, N, K))
    total_energy = exact_F(system) - pristine_system
    # Calculate the energy of formation for each impurity
    impurity_energy =
        map(x -> exact_F(System([x], T, N, K)) - pristine_system, Imps) |> sum
    return (total_energy - impurity_energy)
end

# -TS_I term for the exact diagonalization
function exact_neg_TS_I(system)
    Imps = system.Imps
    T = system.T
    N = system.N
    K = system.K
    pristine_system = exact_neg_TS(System([], T, N, K))
    total_energy = exact_neg_TS(system) - pristine_system
    # Calculate the energy of formation for each impurity
    impurity_energy =
        map(x -> exact_neg_TS(System([x], T, N, K)) - pristine_system, Imps) |>
        sum
    return (total_energy - impurity_energy)
end

## Functions
function Ω(θ, K)
    return √(4 * sin(θ / 2)^2 + K)
end

function Π_jl(z, K, D)

    f_int(θ) = 1 ./ (-1 .+ (Ω(θ, K) ./ z)^2) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π)[1]
    return (res / (2 * π * z^2))
end

function E_I_Integrand(ω, system)
    imps = system.Imps
    K = system.K

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    m = ones(nImps) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    Π = map((x, y) -> Π_jl(ω, K, x.pos - y.pos), imps_mat, imps_mat_T)
    # Matrix to remove couplings between defects
    mask = Diagonal(ones(nImps))
    ΠD = mask .* Π

    Ξ = [
        m_sqrt_inv*Π*m_sqrt_inv ω*m_sqrt_inv*Π*m_sqrt
        ω*m_sqrt*Π*m_sqrt_inv (ω^2)*m_sqrt*Π*m_sqrt+m
    ]
    ΞD = [
        m_sqrt_inv*ΠD*m_sqrt_inv ω*m_sqrt_inv*ΠD*m_sqrt
        ω*m_sqrt*ΠD*m_sqrt_inv (ω^2)*m_sqrt*ΠD*m_sqrt+m
    ]

    unit_mat = ones(2 * length(imps)) |> Diagonal
    Δ_Λ = vcat(map(x -> x.Δ, imps), map(x -> 1 / x.M - 1, imps)) |> Diagonal
    return (
        (det(unit_mat .+ Ξ * Δ_Λ) |> Complex |> log) -
        (det(unit_mat .+ ΞD * Δ_Λ) |> Complex |> log)
    )
end

# E_I
function E_I(system)
    if system.T == 0
        res = quadgk(
            ω -> 2 * real(E_I_Integrand(ω * 1im, system)),
            0,
            Inf,
            maxevals = 1e5,
            rtol = ν,
        )[1]::Float64
        return (res / (2 * π) / 2)
    else
        max_n = Integer(floor(max_omega / (2 * π * system.T)))
        res =
            map(
                n ->
                    system.T *
                    (E_I_Integrand(2 * pi * n * system.T * 1im, system)),
                1:1:1000,
            ) |>
            sum |>
            real
        return res
    end
end

## Zero Matsubara Frequency: -TS
function Π0_jl(K, D)

    f_int(θ) = 1 ./ (Ω(θ, K)^2) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π)[1]
    return (res / (2 * π))
end

# -TS_I
function neg_TS_I(system)
    imps = system.Imps
    K = system.K

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    m = ones(nImps) |> Diagonal
    m_sqrt = sqrt(m)
    m_sqrt_inv = inv(m_sqrt)

    unit_mat = ones(length(imps)) |> Diagonal

    Π = map((x, y) -> Π0_jl(K, x.pos - y.pos), imps_mat, imps_mat_T)
    # Matrix to remove couplings between defects
    mask = Diagonal(ones(nImps))
    ΠD = mask .* Π

    res =
        (
            unit_mat +
            m_sqrt_inv * Π * m_sqrt_inv * Diagonal(map(x -> x.Δ, imps)) |>
            det |>
            log
        ) - (
            unit_mat +
            m_sqrt_inv * ΠD * m_sqrt_inv * Diagonal(map(x -> x.Δ, imps)) |>
            det |>
            log
        )
    return (res * system.T / 2) |> real
end

Imp_1 = Impurity(1, 2, 0)
Imp_2 = Impurity(3, 2, 0)

@time exact_F_I(System([Imp_1, Imp_2], 0.02, 300, 1))
@time E_I(System([Imp_1, Imp_2], 0.02, 300, 1))
# return exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K))
@time E_I_Integrand(10000im, System([Imp_1, Imp_2], 0.4, 300, 0.1)) |> real
