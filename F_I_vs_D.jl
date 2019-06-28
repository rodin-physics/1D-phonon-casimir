using LaTeXStrings
include("Casimir_Effect_Library.jl")

## Parameters

# Separations
Ds = 1 : 1 : 50

# T in units of Ω0
T = 2e-2;

# αs: 1 for infinitely heavy impurity, 0 for the same mass as the chain atoms
αs = [0.25, 0.75, 0.9, 1 - 1e-3]

## Plotting
pgfplots()
plot(
    xaxis = (L"$\ln D$", font(14)),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14)),
    xtickfont = font(12),
    ytickfont = font(12),
    legendfont = font(12),
    legend = :bottomleft,
    )

# Bounding power laws
r_1 = 1 ./ Ds
r_3 =1 ./ (Ds.^3)

for ii = 1 : length(αs)
    α = αs[ii];
    E0 = Energy_Integral(1, α, T)
    r =  map(x -> Energy_Integral(x, α, T) / E0, Ds);
    plot!(log.(Ds), log.(r), linewidth = 2,  color = colors[ii], lab = latexstring(L"$\alpha = $" * string(α)))
end

plot!(log.(Ds), log.(r_1), linewidth = 2,color = colors[5], line = :dash, lab = :false)
plot!(log.(Ds), log.(r_3), linewidth = 2,color = colors[5], line = :dash, lab = :false)

savefig("F_I_T2e_2.pdf")
