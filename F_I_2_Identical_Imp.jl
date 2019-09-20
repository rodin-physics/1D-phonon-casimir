include("Casimir_Effect_Library.jl")

## Parameters
# Separations
Ds = 1 : 1 : 50         # Distances used for the analytic calculation
Ds_Exact = 2 : 3 : 50   # Distances used for the exact diagonalization
N = 1000                # Chain length for the exact diagonalization
# T in units of Ω0
T = 1e-12;

# Impurity masses in units of m (chain masses)
Ms = [1/3, 4/3, 4, 10, 1000]

## Plotting
pyplot()
plot(
    xaxis = (L"$\ln\, D$", font(14, "Serif")),
    yaxis = (L"$\ln \left(F_I/F_I^1\right) $", font(14, "Serif")),
    xtickfont = font(12, "Serif"),
    ytickfont = font(12, "Serif"),
    legendfont = font(12, "Serif"),
    legend = :bottomleft,
    yticks = -10 : 5 : 0
    )

# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds.^3)

for ii = 1 : length(Ms)
    println(ii)
    M = Ms[ii];
    f(x) = F_I([Impurity(0, M), Impurity(x, M)], T)
    # F_I for adjacent impurities
    E0 = f(1)
    α = 1 - 1 / M
    # F_I divided by F_I at D = 1
    r =  map(f , Ds) ./ E0;
    plot!(log.(Ds), log.(r),
        linewidth = 2,
        color = colors[ii],
        lab = latexstring(L"$\alpha = $" * string(α))
        )
end

for ii = 1 : length(Ms)
    println(ii)
    M = Ms[ii];
    # Energy for maximally separated impurities in a finite-length chain
    E_halfway = Exact_Energy(N, M, M, floor(Int, N / 2), T)
    # F_I for adjacent impurities
    E0 = Exact_Energy(N, M, M, 1, T) - E_halfway
    # F_I divided by F_I at D = 1
    r =  map(x -> Exact_Energy(N, M, M, x, T) - E_halfway, Ds_Exact) ./ E0;
    plot!(log.(Ds_Exact), log.(r),
        color = colors[ii],
        lab = "",
        markershape = :circle
        )
end

plot!(
    log.(Ds), log.(r_1),
    linewidth = 2,
    color = RGB(50/255,50/255,50/255),
    line = :dash,
    lab = ""
    )

plot!(
    log.(Ds), log.(r_3),
    linewidth = 2,
    color = RGB(50/255,50/255,50/255),
    line = :dash,
    lab = ""
    )

savefig("F_I_T1e_12.pdf")
