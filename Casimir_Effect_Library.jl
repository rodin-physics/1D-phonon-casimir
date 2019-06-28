using QuadGK
using Plots

## Parameters
ν = 1e-5;       # Relative tolerance for integration

## Functions
# Bose-Einstein distribution
function nB(x)
    return (1 ./ (exp(x) - 1))
end

# Kernel used in the integration to compute F_I for two impurities
function Energy_Kernel(θ, D, α)
    num = exp(4im * D * θ);
    den = 1 + 1im * (1 / α - 1) * cot(θ);
    return (log( 1 - num / den^2 ))
end

# F_I for two impurirties. The temperature parameter t = T / Ω
function Energy_Integral(D, α, t)
    int_res = quadgk(θ -> cos(θ) * Energy_Kernel(θ, D, α) .* (0.5 + nB(sin(θ)/t)), 0, π / 2, rtol = ν)
    return imag(int_res[1] / π)
end

# Colors for plotting
colors = [RGB(215/255,67/255,84/255)
        , RGB(106/255,178/255,71/255)
        , RGB(100/255,101/255,218/255)
        , RGB(169/255,89/255,201/255)
        , RGB(209/255,135/255,46/255)
         ]
