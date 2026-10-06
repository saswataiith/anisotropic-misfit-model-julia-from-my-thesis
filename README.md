# Anisotropic misfit model from my thesis — Julia

This repository contains a Julia implementation of my hexagonal-to-orthorhombic phase-field model. I developed the original codes for my thesis, *Evolution of Multivariant Microstructures with Anisotropic Misfit: A Phase Field Study*. I also used these codes to study symmetry breaking and shape transitions.

The code uses the CPU. It runs on macOS, Linux and Windows with Julia and FFTW; a GPU is not required. GPU acceleration is not implemented.

## Run the example

Install Julia, open a terminal in this repository, and run:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. hex_to_orth.jl
julia --project=. checks.jl
```

The example uses 32 × 32 points, spacing 0.5, timestep 0.001 and 100 steps. The initial mean composition is 0.3. Small Gaussian fluctuations are generated with seed 494. This is a short numerical check, not a simulation of a developed precipitate or a reproduction of a thesis figure. Change the initial fields and run parameters for a research calculation.

The `results` directory contains composition and three variant fields as CSV arrays. Rows follow x and columns follow y. `history.csv` contains time, mean composition and total free energy. Parameters are collected at the beginning of `hex_to_orth.jl`; the functions then construct the grid, calculate elastic kernels, calculate derivatives, advance the fields and save the results.

## Model

Composition c is conserved. The three fields η₁, η₂ and η₃ describe the structural variants. All values are dimensionless: the recovered parameter file does not provide a conversion to SI length, time or energy. Strain and composition are dimensionless as well.

Let qₐ=ηₐ². The local free-energy density used in the selected historical routine `evolve12.c` is

f=A₁c²/2+A₂(1−c)∑qₐ+A₄₁∑qₐ²+A₄₂∑ₐ<ᵦqₐqᵦ+A₆₁∑qₐ³+A₆₂q₁q₂q₃.

A₁ and A₂ set the quadratic terms and composition–variant coupling. A₄₁ and A₄₂ set the fourth-order terms; A₆₁ and A₆₂ set the sixth-order terms. The gradient energy is κc|∇c|²+κη∑|∇ηₐ|². Each κ is its corresponding gradient-energy coefficient.

The transformation strain is ε⁰=εc c+∑εₐηₐ². Here εc is the composition strain tensor; εₐ is the strain tensor for variant a. `eps_c` gives the in-plane dilatational composition strain. `eps_eta` gives the first principal variant strain; `tetra` gives the ratio of the second principal strain to the first. The other two variants are rotations by ±60°.

The elastic moduli are homogeneous. `mu`, `nu` and `anisotropy` use the same conversion to cubic stiffnesses as my C routine `calc_bn.c`. These parameters specify the stiffness tensor C. For each nonzero Fourier direction n, the acoustic tensor is Qᵢₗ=Cᵢⱼₖₗnⱼnₖ. With σₐ=C:εₐ and traction aₐ=σₐn, the interaction kernel is Bₐᵦ=εₐ:C:εᵦ−aₐ·Q⁻¹aᵦ. The code uses all four strain contributions: composition and three variants.

The grid is two-dimensional, with n_z=0 and ε⁰_zz=0. The stiffness calculation retains the three-dimensional tensor, matching the historical construction. The zero Fourier mode is omitted, corresponding to relaxed homogeneous strain. A prescribed mean strain would require its additional homogeneous energy term.

The transport equations are ∂c/∂t=M∇²(δF/δc) and ∂ηₐ/∂t=−LδF/δηₐ. M is composition mobility; L is variant relaxation coefficient. Periodic boundaries and semi-implicit Fourier updates are used. Gradient terms are implicit and local/elastic terms are explicit. This does not establish unconditional energy stability.

## What differs from the historical program

This implementation selects the `evolve12.c` equations; it does not combine all historical variants. It uses equal variant gradient coefficients and relaxation coefficients, a constant timestep, a simple reproducible initial field and no sustained noise. Julia's random sequence differs from Numerical Recipes and GSL.

The elastic kernels are explicitly symmetrized between conjugate Fourier modes on even-grid Nyquist lines, so real fields have a real convolution and a consistent discrete energy. This is a documented numerical change from assigning a single sign on those lines in the historical code. No claim of bitwise agreement with the old program is made.

`check-results.txt` records the completed checks. They cover energy derivatives, nonnegative elastic quadratic forms, finite fields, conservation and decreasing energy in the short example. These are implementation checks; grid convergence and published shape-transition results remain to be reproduced.

## Particle moments and shape parameter

My original utility is retained in `reference/calc_singleppt.c`. It refers to Numerical Recipes headers, but those third-party implementations are not included in this public repository. The new Julia utility needs only Julia's standard library.

```julia
include("utilities/particle_moments.jl")
using DelimitedFiles
c = readdlm("results/composition.csv", ',')
eta = readdlm("results/variant-1.csv", ',')
# Use a developed particle field; the short noise example has no selected particle.
result = particle_moments(c, eta; dx=0.5, dy=0.5)
println(result.shape_parameter)
```

The selected region satisfies c≥0.5 and η≥0.5. Use one isolated precipitate. The utility calculates its area, area-centroid, equivalent circular radius, second-moment tensor, principal moments and long-axis angle. It also reads the original C binary files through `read_c_field(path,nx,ny)`.

Coordinates x and y are distances from the area-centroid. For uniform density, Ixx=∫y²dA, Iyy=∫x²dA and Ixy=−∫xy dA. These are area moments, with units length⁴ if the spacing has physical length units. Multiplying by mass per unit area would give mass moments of inertia. Area has units length²; the equivalent radius has units length.

The shape parameter is (Imax−Imin)/(Imax+Imin), where Imax and Imin are the principal moments. It is dimensionless. A circle gives zero. A uniform ellipse with semiaxes a and b gives |a²−b²|/(a²+b²). Unequal moments identify elongation. Equal moments do not prove that the outline is circular: some noncircular shapes also have equal moments. Use contours or higher-order measures to examine those changes.

The angle is in degrees, measured from +x modulo 180°. It is undefined when the moments are equal, so the utility returns NaN. The mask must not cross a periodic boundary; shift the isolated particle into the box first. It does not automatically separate several particles.

The Julia version deliberately corrects several assumptions in the original utility: it measures about the particle centroid instead of the box centre, applies dx and dy to distances, and obtains the orientation from the symmetric tensor rather than a one-argument arctangent. Circle, ellipse, translation and spacing-scaling checks are included.

## Historical source

The original C archive and the GSL working version are preserved in [my historical repository](https://github.com/saswataiith/anisotropic-misfit-model-from-my-thesis), which remains private because it contains third-party Numerical Recipes source. No new licence has been imposed on that archive. The new Julia implementation does not depend on Numerical Recipes.

The `gsl/` directory contains the previously tested C working version, with Numerical Recipes replaced by GSL. Build with `make -C gsl/src_serial` after installing FFTW and GSL (on macOS: `brew install fftw gsl`). It retains the historical input format and three separately specified variant mobilities and gradient coefficients. The Julia version is the simpler implementation described above. Selected original model routines are included in `reference/` for comparison.
