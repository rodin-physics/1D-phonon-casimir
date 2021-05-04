N = 3000;
ΩT = 1e1;
K = 0.5;
s = System([], ΩT, N, K)
eig = exact_modes(s)

nPts = 40000;
numerical_displacement = map(x -> displacement([1], eig, ΩT, [1])[1], 1:nPts)


# Analytic standard deviation
analytic_std = pristine_correlation(K, 0, ΩT) |> sqrt
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
