using CairoMakie
using DelimitedFiles
using LinearAlgebra
using ProgressMeter
using QuadGK

## Parameters
const ν = 1e-5;             # Relative tolerance for integration
const η = 1e-5;             # Small number used for i0
const α = 1e-9;             # Small number for absolute tolerance
const max_omega = 100;      # Large number used to truncate the finite
# temperature Matsubara sum
const NumEvals = 1e7;       # Maximum number of evaluations in quadgk
