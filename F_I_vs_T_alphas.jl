using LaTeXStrings
include("Casimir_Effect_Library.jl")

## Parameters
nPts = 200;

# T in units of Ω0
T_min = 1e-12;
T_max = 2e-1;
Ts = range(T_min, stop = T_max, length = nPts)

αs = [0.25, 0.75, 0.9, 1 - 1e-3]
D = 1;

## Plotting

pgfplots()
plot(
    xaxis = (L"$T/\Omega_0$", font(14)),
    yaxis = (L"$ \left|F_I^1/\Omega^0\right|$", font(14)),
    xtickfont = font(12),
    ytickfont = font(12),
    legendfont = font(12),
    legend = :topright,
    )

for ii = 1 : length(αs)
    α = αs[ii];
    r =  map(x -> Energy_Integral(D, α, x), Ts);
    plot!(Ts, abs.(r),
        linewidth = 2,
        color = colors[ii],
        lab = latexstring(L"$\alpha = $" * string(α))
        )
end

savefig("Finite_T_Energy_alphas.pdf")
