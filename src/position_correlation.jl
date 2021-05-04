include("general.jl")

function pristine_correlation(T, D, K)

    f_int(θ) = 1 ./ Ω(θ, K) .* (nB(Ω(θ, K) / T) + 1 / 2) * exp(1im * θ * D)

    res = quadgk(f_int, 0, 2 * π, atol = 1e-9)[1]
    return (res / (2 * π))
end

function correlation_correction_Integrand(ω, system, j, l)
    imps = system.Imps
    K = system.K

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    m = ones(nImps) |> Diagonal

    Π = map((x, y) -> Π_jl(ω, K, x.pos - y.pos), imps_mat, imps_mat_T)

    Ξ = [
        Π ω*Π
        ω*Π (ω^2)*Π+m
    ]

    unit_mat = ones(2 * length(imps)) |> Diagonal
    Δ_Λ = vcat(map(x -> x.Δ, imps), map(x -> 1 / x.M - 1, imps)) |> Diagonal
    prop = Δ_Λ * inv(unit_mat .+ Ξ * Δ_Λ)

    # Calculate the coupling between the mass of interest and the defects
    Π_j = map(x -> Π_jl(ω, K, x.pos - j), imps)
    Π_l = map(x -> Π_jl(ω, K, x.pos - l), imps)

    left_prop = [Π_j; ω * Π_j] |> permutedims
    right_prop = [Π_l; ω * Π_l]

    return (left_prop*prop*right_prop)[1]
end


function correlation_correction_zero_Matsubara(system, j, l)
    imps = system.Imps
    K = system.K

    nImps = length(imps)
    imps_mat = repeat(imps, 1, nImps)
    imps_mat_T = permutedims(imps_mat)

    m = ones(nImps) |> Diagonal

    Π = map((x, y) -> Π0_jl(K, x.pos - y.pos), imps_mat, imps_mat_T)

    Ξ = [
        Π zeros(nImps, nImps)
        zeros(nImps, nImps) m
    ]

    unit_mat = (ones(2 * length(imps)) |> Diagonal)
    Δ_Λ = vcat(map(x -> x.Δ, imps), map(x -> 1 / x.M - 1, imps)) |> diagm

    prop = Δ_Λ * inv(unit_mat .+ Ξ * Δ_Λ)


    # Calculate the coupling between the mass of interest and the defects
    Π_j = map(x -> Π0_jl(K, x.pos - j), imps)
    Π_l = map(x -> Π0_jl(K, x.pos - l), imps)

    # left_prop = [Π_j; 0 * Π_j] |> permutedims
    # right_prop = [Π_l; 0 * Π_l]
    left_prop = [Π_j; zeros(nImps)] |> permutedims
    right_prop = [Π_l; zeros(nImps)]

    return (left_prop*prop*right_prop)[1]
end

function correlation_correction(system, j, l)
    if system.T == 0
        res = quadgk(
            ω ->
                -2 * real(
                    correlation_correction_Integrand(ω * 1im, system, j, l),
                ),
            0,
            Inf,
            maxevals = NumEvals,
            rtol = ν,
            atol = 1e-5,
        )[1]::Float64
        return (res / (2 * π))
    else
        max_n = Integer(floor(max_omega / (2 * π * system.T)))
        res =
            map(
                n ->
                    -2 *
                    system.T *
                    (correlation_correction_Integrand(
                        2 * pi * n * system.T * 1im,
                        system,
                        j,
                        l,
                    )),
                1:1:1000,
            ) |>
            sum |>
            real
        return (
            res -
            correlation_correction_zero_Matsubara(system, j, l) * system.T
        )
    end
end
