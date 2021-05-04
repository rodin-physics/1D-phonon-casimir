include("../src/QFT_energy.jl")
include("../src/exact_system.jl")
include("../src/position_correlation.jl")

# Energy testing
N = 2000;
K = 0.5;
T = 1e-1;
Imp_1 = Impurity(1, 5, 0)
Imp_2 = Impurity(2, 8, 0)
Imp_3 = Impurity(4, 4, 0)

s = System([Imp_1, Imp_2, Imp_3], T, N, K)

exact_F_I(s) - exact_neg_TS_I(s)
E_I(s)

# Displacement testing
N = 2000;
K = 0.5;
T = 0.1;
loc = 4;

Imp_1 = Impurity(1, 5, 0)
Imp_2 = Impurity(2, 1, 100)
Imp_3 = Impurity(4, 2, 2)

s = System([Imp_1, Imp_2, Imp_3], T, N, K)
eig = exact_modes(s)
nPts = 40000;
numerical_displacement = map(x -> displacement([loc], eig, s)[1], 1:nPts)

analytic_std =
    pristine_correlation(T, 0, K) + correlation_correction(s, loc, loc) |>
    sqrt |>
    real
disp_points = range(-5, 5, length = 100);
analytic_distribution =
    exp.(-disp_points .^ 2 / 2 / analytic_std^2) / (analytic_std * √(2 * π))

fig = Figure(resolution = (1800, 1800))
ax =
    fig[1, 1] = Axis(
        fig,
        xlabel = "δ",
        ylabel = "P(δ)",
        xlabelpadding = 0,
        ylabelpadding = 0,
        xlabelsize = 12,
        ylabelsize = 12,
        xticklabelsize = 12,
        yticklabelsize = 12,
        aspect = AxisAspect(1),
        xticklabelfont = "serif-roman",
        yticklabelfont = "serif-roman",
        xlabelfont = "serif-italic",
        ylabelfont = "serif-italic",
    )
lines!(disp_points, analytic_distribution, linewidth = 2)

density!(
    numerical_displacement,
    npoints = 100,
    color = (:red, 0.3),
    strokecolor = :red,
)
fig
