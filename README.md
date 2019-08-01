# 1D-phonon-casimir

This repository includes the code used to calculate and plot the results in...

The file `Casimir_Effect_Library.jl` contains the functions used in other files:
* `F_I_vs_D.jl`is used to compute the interaction energy between two impurities at a given $T$ for several values of $\alpha$
* `F_I_vs_T_Ds.jl` calculates the interaction energy as a function of $T$ at a fixed $\alpha$ for several separations
* `F_I_vs_T_alphas.jl` computes the interaction energy as a function of $T$ at a fixed $D$ for several $\alpha$s
* `Energy_Correction.jl` is used to calculate the energy correction for each phonon momentum $\theta$ at a fixed $D$ for several values of $\alpha$
