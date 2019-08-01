using LaTeXStrings
include("Casimir_Effect_Library.jl")

## Parameters

nPts = 3000;
θs = range(0, π / 2, length = nPts)

# αs: 1 for infinitely heavy impurity, 0 for the same mass as the chain atoms
αs = [0.25, 0.75, 0.9, 0.99]
D = 1;

## Plotting

pgfplots()
plot(
    xaxis = (L"$\theta$", font(14)),
    yaxis = (L"$\cos\theta\times\mathrm{Im}\left[\dots\right]/\pi $", font(14)),
    xtickfont = font(12),
    ytickfont = font(12),
    legendfont = font(12),
    legend = :bottomright,
    ylims = (-0.5, 0.5)
    )


for ii = 1 : length(αs)
    α = αs[ii];
    res = imag.(Energy_Kernel.(θs, D, α)) ./ π .* cos.(θs)
    plot!(θs, res,
        linewidth = 2,
        color = colors[ii],
        lab = latexstring(L"$\alpha = $" * string(α))
        )
end

# Plot the position of the isolated dimer mode as a vertical dashed line
for ii = 1 : length(αs)
    α = αs[ii];
    xs = asin(sqrt((1 - α) / 2)) .* [1, 1]
    ys = [-1, 1]
    plot!(xs, ys,
        linewidth = 2,
        color = colors[ii],
        line = :dash,
        lab = :false
        )
end

savefig("Dimer_Mode_D1.pdf")
