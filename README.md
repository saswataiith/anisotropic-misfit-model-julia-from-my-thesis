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

The example uses 32 × 32 points, spacing 0.5, timestep 0.001 and 100 steps. The initial mean composition is 0.3. Small Gaussian fluctuations are generated with seed 494. I retain sustained noise: the default amplitudes are 0.1 for c and each variant, with the original event schedule. This is a short numerical check, not a simulation of a developed precipitate or a reproduction of a thesis figure. Change the initial fields and run parameters for a research calculation.

The `results` directory contains composition and three variant fields as CSV arrays. Rows follow x and columns follow y. `history.csv` contains time, mean composition, total free energy after noise, a noise-event flag, and energy before noise. The initial row has no pre-noise energy (NaN). Parameters are collected at the beginning of `hex_to_orth.jl`; the functions then construct the grid, calculate elastic kernels, calculate derivatives, advance the fields and save the results.

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

This implementation selects the `evolve12.c` equations; it does not combine all historical variants. It supports separately specified variant gradient coefficients and relaxation coefficients, two timesteps with a switch count, and the sustained-noise block described below. The initial random composition uses the original uniform distribution with its mean removed; each variant starts with Gaussian noise. Julia's random sequence differs from Numerical Recipes and GSL.

The elastic kernels are explicitly symmetrized between conjugate Fourier modes on even-grid Nyquist lines, so real fields have a real convolution and a consistent discrete energy. This is a documented numerical change from assigning a single sign on those lines in the historical code. No claim of bitwise agreement with the old program is made.

`check-results.txt` records the completed checks. They cover energy derivatives, nonnegative elastic quadratic forms, finite fields, conservation and decreasing energy in a deterministic short run with noise amplitudes explicitly set to zero. Separate tests exercise the noise schedule, conserved composition perturbations, all three nonconserved variant perturbations, seeded repeatability, a noisy evolution run, timestep switching and separate variant coefficients. These are implementation checks; grid convergence and published shape-transition results remain to be reproduced.

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

The recovered model sources and GSL working version are publicly available in [my historical C repository](https://github.com/saswataiith/hex-to-orth-historical-codes). The complete byte-preserved archive is in [my private archive repository](https://github.com/saswataiith/anisotropic-misfit-model-from-my-thesis), which remains private because it contains third-party Numerical Recipes source. No new licence has been imposed on that archive. The new Julia implementation does not depend on Numerical Recipes.

The `gsl/` directory contains the previously tested C working version, with Numerical Recipes replaced by GSL. Build with `make -C gsl/src_serial` after installing FFTW and GSL (on macOS: `brew install fftw gsl`). It retains the historical input format and three separately specified variant mobilities and gradient coefficients. The Julia version is the simpler implementation described above. Selected original model routines are included in `reference/` for comparison.

## Sustained noise: retained explicitly

I follow the noise block in `evolve12.c`. After each deterministic update, I add a Gaussian perturbation to c and subtract its spatial mean. This preserves the total composition. I independently add Gaussian perturbations to eta1, eta2 and eta3 without subtracting their means, because these fields are nonconserved.

The event condition is `count <= initial_steps || count % interval == 0`. Count starts at zero, as in the C code. The default interval is 2500. Therefore even initial_steps=0 permits the first event at count zero. Setting initial_steps=0 alone does not turn off noise.

The amplitudes are increments per noise event, in the same dimensionless field units as c and eta. I do not multiply them by sqrt(dt) or reinterpret them as a thermally calibrated fluctuation–dissipation model. Adding noise can increase energy: decreasing energy is checked separately with all amplitudes zero. I do not clip noisy fields.

```julia
include("hex_to_orth.jl")
noise = Noise(amplitude_c=0.02, amplitude_eta=(0.01,0.02,0.03),
              initial_steps=100, interval=2500)
run_example(noise=noise)

# I explicitly disable noise only for a deterministic comparison.
quiet = Noise(amplitude_c=0.0, amplitude_eta=(0.0,0.0,0.0))
run_example(noise=quiet, directory="results-deterministic")
```

`run_simulation` accepts the complete composition and three variant arrays, separate parameters, a seeded random generator, dx and dy, dt1 and dt2, time_to_change, initial_count and initial_time. It switches to dt2 when count > time_to_change, as in the selected C routine. I make exactly the requested number of updates and do not reproduce the old loop's extra final update.

To run the supplied historical small parameter file:

```julia
include("hex_to_orth.jl")
run_from_input("gsl/examples/smoke32/InputParams")
```

The parameter reader supports the random initialization branch (flag=0, initcount=0). For other historical initializers, use the complete C versions in the historical repository, or load their saved fields with `read_c_field` and pass those arrays to `run_simulation`. The Julia generator reproduces its own seeds, not the old Numerical Recipes or GSL random sequence. For an exact continuation of a noisy Julia trajectory, retain the random generator state as well as the fields.

## What I tested and what I have not tested

- Energy variation agrees with analytical local, gradient and elastic derivatives.
- Elastic kernels pass a positive-semidefinite check and a known isotropic dilatational reference check.
- A 32² deterministic run of 100 steps preserves mean composition within 2.3e−16 and decreases energy.
- Noise events follow the original count schedule. Composition noise has zero spatial mean; each variant receives a separate nonconserved perturbation. Identical seeds repeat the same noise arrays.
- A 16² run with repeated noise events preserves mean composition within 1.2e−16.
- Timestep switching and unequal variant gradient/relaxation coefficients are exercised.
- The shape utility passes circle, ellipse, translation, spacing-scaling and empty-region checks.

These are functional and equation-consistency checks. I have not reproduced a published symmetry-breaking shape transition, established grid convergence, or validated every historical evolution routine and initializer. The complete recovered model sources are retained publicly; I do not claim every historical branch has been translated into Julia.

## Development of the auxiliary-variable models

[Read my model history](MODEL-HISTORY.md), including cA in my PRL ternary model, phi for boundary conditions in my PCCP model, and the free-energy interpolation explored in Sandeep's thesis, Chapter 3 Eq. (3.77) and Chapter 6 Section 6.1.

## Using my code

My original source is available under MIT. See [licence scope and dependency terms](LICENSING.md) and the [licence](LICENSE).
