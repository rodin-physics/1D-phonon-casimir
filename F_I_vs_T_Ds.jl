using LaTeXStrings
include("Casimir_Effect_Library.jl")

## Parameters
nPts = 200;

# T in units of Ω0
T_min = 1e-12;
T_max = 1e-1;
Ts = range(T_min, stop = T_max, length = nPts)

Ds = [1, 5, 10, 20];
α = 0.5;

# Vertical line at T = 0.02
xs = [0.02, 0.02]
ys = [0, 1]

## Plotting
pgfplots()
plot(
    xaxis = (L"$T/\Omega_0$", font(14)),
    yaxis = (L"$ F_I/F_I^{T=0}$", font(14)),
    xtickfont = font(12),
    ytickfont = font(12),
    legendfont = font(12),
    legend = :bottomright,
    )

for ii = 1 : length(Ds)
    D = Ds[ii];
    # F_I for adjacent impurities
    E0 = Energy_Integral(1, α, T)
    # F_I divided by F_I at D = 1
    r =  map(x -> Energy_Integral(D, α, x) / E0, Ts);
    plot!(Ts, r,
        linewidth = 2,
        color = colors[ii],
        lab = latexstring(L"$D = $" * string(D))
        )
end

# Add the vertical line
plot!(xs, ys,
    linewidth = 2,
    color = colors[5],
    line = :dash,
    lab = :false
    )

savefig("Finite_T_Energy_Ds.pdf")
