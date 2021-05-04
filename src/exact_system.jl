include("general.jl")

# Function to obtain the modes of the system with defects
function exact_modes(system)
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
    M_mat_inv_sqrt = sqrt.(inv(M_Mat))
    res = eigen(M_mat_inv_sqrt * U_Mat * M_mat_inv_sqrt)
    return res
end

# Exact Diagonalization 1D to obtain the total Free Energy
function exact_F(system)
    T = system.T
    eigs = exact_modes(system)
    Ω2 = eigs.values |> real
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
        Δs[ii.pos] = ii.Δ
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

# Randomly-generated ξq
function ξq(Ωq, ΩT)
    σ = sqrt(coth(Ωq / 2 / ΩT)) / 2
    return σ * (randn() + im * randn())
end

function massof(pos, system)
    Imps = system.Imps
    imp_pos = map(x -> x.pos, Imps)
    imp_mass = map(x -> x.M, Imps)
    mass = imp_mass[imp_pos.==pos]
    if length(mass) == 0
        return 1.0
    else
        return mass[1]
    end
end

function displacement(pos, eigs, T, system)
    Ω2 = eigs.values |> real
    Ωs = sqrt.(abs.(Ω2))
    ξs = map(x -> ξq(x, T), Ωs)
    amp = √(2) .* real(ξs) ./ sqrt.(Ωs)

    res = map(
        x -> amp .* (eigs.vectors)[x, :] ./ sqrt(massof(x, system)) |> sum,
        pos,
    )

    return (res)
end
