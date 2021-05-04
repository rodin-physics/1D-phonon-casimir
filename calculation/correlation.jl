include("../src/position_correlation.jl")
include("../src/exact_system.jl")
# Displacement testing
N = 1500;
K = 1e-3;
T = 0;
loc = 4;

Imp_1 = Impurity(1, 1e12, 0)
Imp_2 = Impurity(5, 1e12, 0)
Imp_3 = Impurity(4, 2, 2)

s = System([Imp_1, Imp_2], T, N, K)
eig = exact_modes(s)

nPts = 200000;
pair_displacement = map(x -> displacement([1, 5], eig, s), 1:nPts)
d = reduce(hcat, pair_displacement)
cor(d[1, :], d[2, :])
pristine_correlation(T, 3, K) + correlation_correction(s, 1, 4)

var(d[1, :])
(pristine_correlation(T, 0, K) + correlation_correction(s, 1, 1))


(pristine_correlation(T, 4, K) + correlation_correction(s, 1, 5)) /
sqrt(pristine_correlation(T, 0, K) + correlation_correction(s, 1, 1)) /
sqrt(pristine_correlation(T, 0, K) + correlation_correction(s, 5, 5))

pristine_correlation(T, 0, K) + correlation_correction(s, 1, 1)
(pristine_correlation(T, 4, K)) / sqrt(pristine_correlation(T, 0, K)) /
sqrt(pristine_correlation(T, 0, K))


(pristine_correlation(50.0, 12, K)) / sqrt(pristine_correlation(50.0, 0, K)) /
sqrt(pristine_correlation(50.0, 0, K))


nPts = 40000;
numerical_displacement = map(x -> displacement([loc], eig, s)[1], 1:nPts)

analytic_std =
    pristine_correlation(T, 0, K) + correlation_correction(s, loc, loc) |>
    sqrt |>
    real
std(numerical_displacement)




pristine_correlation(T, 1, K) + correlation_correction(s, 1, 3)
