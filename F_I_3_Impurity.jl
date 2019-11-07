using LaTeXStrings
include("Casimir_Effect_Library.jl")

## Parameters
# Separations
Ds = 1 : 1 : 10

M_inner = .75;
M_outer = 20;

nTemps = 50;
ln_T_min = log(1e-2);
ln_T_max = log(1e-1);

ln_Ts = range(ln_T_min, ln_T_max, length = nTemps);

ln_Ts_Mat = repeat(ln_Ts, 1, length(Ds))
Ds_Mat = repeat(Ds', nTemps,1)

## Functions
f(r, T) = F_I([Impurity(-r, M_outer), Impurity(0, M_inner), Impurity(r, M_outer)], T)

## Calculation

println("Calculating")
res = map((r, T) ->  f(r, exp(T)), Ds_Mat, ln_Ts_Mat)
res = res * 1e4;
bound = maximum(abs.(res))/10

println("Plotting")
## Plotting
pyplot()

heatmap(Ds, ln_Ts, res,
        # leg = false,
        # aspect_ratio = 1,
        # size = (500,400),
        xaxis = (L"D"),
        yaxis = (L"\ln \left(T / \Omega\right)"),
        xtickfont = font(12, "Serif"),
        ytickfont = font(12, "Serif"),
        color = :coolwarm,
        clim = (-bound, bound),
        colorbar = true,
        colorbar_title = L"$F_I/\Omega\times 10^4$",
        # legendfont = "Serif"
        )

println("Saving")
savefig("F_I_3_Imp.pdf")
