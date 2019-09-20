include("Casimir_Effect_Library.jl")

## Parameters
# Separations
Ds = 1 : 1 : 50        # Distances used for the analytic calculation

# Impurity masses in units of m (chain masses)
M = 15.0

cluster_size = [2, 5, 10, 20]

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

for ii = 1 : length(cluster_size)
    println(ii)
    cl = cluster_size[ii]
    E_cluster = F_I_T0(impurity_cluster(cl, 0, M))
    println(E_cluster)
    f(x) = F_I_cluster(x, cl, M) - 2 * E_cluster
    # Cluster Energy
    # F_I for adjacent clusters
    E0 = f(1)
    # F_I divided by F_I at D = 1
    r =  map(f , Ds) ./ E0;
    plot!(log.(Ds), log.(r),
        linewidth = 2,
        color = colors[ii],
        lab = latexstring(L"$L = $" * string(cl))
        )
end

# Bounding power laws
r_1 = 1 ./ Ds
r_3 = 1 ./ (Ds.^3)

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
savefig("Cluster_Interaction_M15.pdf")
