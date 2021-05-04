include("general.jl")

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
            maxevals = NumEvals,
            rtol = ν,
            atol = α,
        )[1]::Float64
        return (res / (2 * π) / 2)
    else
        max_n = Integer(floor(max_omega / (2 * π * system.T)))
        res =
            map(
                n ->
                    system.T *
                    (E_I_Integrand(2 * pi * n * system.T * 1im, system)),
                1:1:max_n,
            ) |>
            sum |>
            real
        return res
    end
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
