include("exact_diag.jl")

Imp_1 = Impurity(1, 4, 1)
Imp_2 = Impurity(3, 2, 1)

s = System([Imp_1, Imp_2], 0.20, 300, 1)


Imps = s.Imps
imp_pos = map(x -> x.pos, Imps)
imp_mass = map(x -> x.M, Imps)

imp_mass[imp_pos .== 2]

mass = imp_mass[imp_pos .= pos] |> sum

Imps = s.Imps
imp_pos = map(x -> x.pos, Imps)
loc = imp_pos[imp_pos.==2]
@time massof(4, s)

@time exact_F(System([Imp_1, Imp_2], 0.20, 300, 1))
@time exact_F_A(System([Imp_1, Imp_2], 0.20, 300, 1))

@time exact_F_I(System([Imp_1, Imp_2], 0.20, 300, 1)) -
      exact_neg_TS_I(System([Imp_1, Imp_2], 0.20, 300, 1))
@time E_I(System([Imp_1, Imp_2], 0.20, 300, 1))
# return exact_F_I(System(Ms, [Imp_1, Imp_2], 0, N, K))
@time E_I_Integrand(10000im, System([Imp_1, Imp_2], 0.4, 300, 0.1)) |> real

@time Π_jl(10im, 1, 1)
@time ΠΩ_jl(10im, 1, 1)
@time Π_Ω_jl(10im, 1, 1)

Ω(π / 2, 1)
c = Π_Ω_jl(10im, 1, 1)

c |> real




Imp_1 = Impurity(1, 1, 1)
Imp_2 = Impurity(5, 1, 100)
pristine_correlation(.0001, 1, 0)
@time correlation_correction(System([Imp_1], 1, 300, 1), 1, 1)



Imp_2 = Impurity(2, 2, 1)


pristine_correlation(1, 0, 0.10)

map(
    x -> correlation_correction_Integrand(
        2 * pi * 0.1 * x * im,
        System([Imp_1, Imp_2], 0.1, 300, 1),
        1,
        1,
    ),
    1:30,
)
@time correlation_correction_Integrand(-1im, System([Imp_1], 0, 300, 1), 1, 1)


1 / 9000.200 - 1
