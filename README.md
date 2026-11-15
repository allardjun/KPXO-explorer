# KP-X-O Explorer

Julia scripts for exploring a kinetic proofreading model of T cell activation.
The model predicts the fraction of T cells that turn on CD69, a marker of activation, from two properties of the peptide they see: its concentration and its binding strength.

The main figure is a 3d plot of that fraction against pMHC concentration, with one dose-response curve for each of five peptides, placed by their dissociation constant KD.
The floor of the plot shows each peptide's potency, and how steeply potency falls as binding weakens.
An interactive version adds a slider for each of the model's eleven parameters.

A browser version of the same model, which needs no installation, runs at <https://allardjun.github.io/261023/>.

## The model

The model is an extension of Lever et al. 2014 ([Nature Reviews Immunology](https://doi.org/10.1038/nri3728)).
Every step is a formula, so no equations have to be solved numerically:

1. **Binding.**
   T cell receptors and peptide–MHC (pMHC) on the opposing cell bind and unbind.
   The density of bound complexes is the smaller root of a quadratic.
2. **Kinetic proofreading.**
   A bound complex must pass `nKP` signaling steps before it unbinds.
   Weak binders unbind sooner, so fewer of their complexes finish.
3. **Switch (X).**
   A Hill function of the fully proofread complexes turns the signal on inside the cell.
4. **Output (O1).**
   A saturating function of the switch gives the fraction of cells that are CD69+.

Units are µM for pMHC concentration and KD, µm⁻² for surface densities, and s⁻¹ for rates.
The parameter `g1` in the code is g̃1 = g1·kon3d, so the off-rate is `g1·KD`.

Seven of the parameters (`fA`, `alpha`, `TT`, `kon2d`, `g1`, `kP`, `lambdaX`) enter the output only through four combinations, which the slider figure displays:

| Combination | Formula |
|---|---|
| x* | TT / α |
| y* | λX / (fA · TT) |
| KD Lauf | TT · kon2d / g̃1 |
| KD McK | kP / g̃1 |

So many different slider settings give the same curves (i.e., there is degeneracy in the model).

### Potency and discrimination power

Potency P15 is the pMHC concentration at which 15% of cells are CD69+.
The discrimination power is the slope of a straight-line fit of log10 P15 against log10 KD over the peptide panel, as defined by Pettmann et al. 2021 ([eLife](https://doi.org/10.7554/eLife.67092)).
Peptides that never reach 15% are left out of the fit.

### Parameter sets

[src/model.jl](src/model.jl) holds two parameter sets, and the slider figure has a button for each.

- `PARAMS_FITTED` comes from one fit of the model to CD69 data.
  The fitting method is not yet published, so treat these values as a starting point, not as a result.
  All the scripts start from this set.
- `PARAMS_NOMINAL` is the initial guess for that fit.

## Installing

You need Julia 1.13 or later.
[juliaup](https://github.com/JuliaLang/juliaup) is the easiest way to get it.

From the repository root, install the exact package versions recorded in `Manifest.toml`:

```sh
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

The first run takes several minutes, since the Makie plotting packages are large and Julia precompiles them.

## Running

Run each script from the repository root.

| Script | What it does | Plotting package |
|---|---|---|
| `scripts/save_cd69_3d.jl` | Writes the figure to `figures/cd69_3d.pdf` | CairoMakie |
| `scripts/interactive_cd69_3d.jl` | Opens the figure in a window you can rotate | GLMakie |
| `scripts/interactive_cd69_3d_web.jl` | Opens the figure in a browser tab you can rotate | WGLMakie |
| `scripts/interactive_cd69_sliders.jl` | Opens the figure with a slider for each parameter, in a browser tab | WGLMakie |

For example:

```sh
julia --project=. scripts/interactive_cd69_sliders.jl
```

In the interactive views, left-drag rotates, scroll zooms, and right-drag moves the view.
From a terminal, the GLMakie script waits until you close its window, and the browser scripts wait until you press Enter.

You can also `include` a script from a Julia REPL started with `--project=.`, for example the one in VS Code.
The script then returns at once and leaves `fig` and `ax` in the REPL, plus `params` for the slider figure.
`params[]` is the parameter set the sliders show now, and the "Print to REPL" button prints it.

The sliders work only while Julia is running, since Julia computes the curves.

### Choosing a plotting package

- **CairoMakie** draws static files and needs no display, so it works on a remote or headless machine.
- **WGLMakie** draws in a web browser.
  It is the most dependable interactive option.
- **GLMakie** opens a native window and needs OpenGL and a display.
  On some macOS setups, for example with an external monitor, it lags while rotating; use the WGLMakie scripts there.

## Layout

```
Project.toml, Manifest.toml   the Julia environment
src/model.jl                  the model, potency and discrimination power, with no dependencies
src/plot_cd69_3d.jl           the 3d figure, for any Makie backend
src/cd69_sliders.jl           the figure with sliders
scripts/                      the scripts to run (see above)
figures/                      generated figures, ignored by git
```

The files in `src/` only define functions and constants; the scripts in `scripts/` load them and do the work.
`src/model.jl` needs no packages, so you can `include` it on its own to compute with the model.

## Sources

- Model: an extension of Lever et al., Nature Reviews Immunology 2014, [doi:10.1038/nri3728](https://doi.org/10.1038/nri3728).
- Discrimination power: Pettmann et al., eLife 2021, [doi:10.7554/eLife.67092](https://doi.org/10.7554/eLife.67092).
- Peptide affinities: Aleksic et al., Immunity 2010.
