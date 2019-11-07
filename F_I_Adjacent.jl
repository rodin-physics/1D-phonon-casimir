include("Casimir_Effect_Library.jl")

## Parameters
nPts = 201;
α_min = -5;
α_max = .999;
αs = range(α_min, α_max, length = nPts)
T = 1e-12;

αs_Mat = repeat(αs, 1, nPts)
αs_Mat_T = permutedims(αs_Mat)

bound = 2e-4;

# Compute the result
res = map((x, y) -> F_I([Impurity(0, 1 / (1 - x)), Impurity(D, 1 / (1 - y))], T), αs_Mat, αs_Mat_T)

res = res .* 1e4
bound = bound * 1e4

## Plotting
pyplot()

heatmap(αs, αs, res,
        leg = false,
        aspect_ratio = 1,
        size = (500,400),
        xaxis = (L"\alpha_1"),
        yaxis = (L"\alpha_2"),
        xtickfont = font(12, "Serif"),
        ytickfont = font(12, "Serif"),
        color = :coolwarm,
        clim = (-bound, bound),
        colorbar = true,
        colorbar_title = L"$F_I^1/\Omega\times 10^4$",
        legendfont = "Serif"
        )

savefig("F_I_2_Alpha.pdf")
