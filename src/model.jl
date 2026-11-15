# Base KP-X-O model of T cell activation (closed form, no equations to solve).
# Units: pMHC concentration x in µM, KD in µM, surface densities in µm⁻², rates in s⁻¹.

# Parameter sets are NamedTuples with the same names as the FlexiModel code.
# In the code, `g1` is the manuscript's g̃1 = g1·kon3d, so that koff = g1·KD.

# Nominal values from `ig_params` in define_model_base.jl (the initial guess for learning).
const PARAMS_NOMINAL = (
    fA = 1.0,
    alpha = 1.0,    # µm⁻² per µM
    TT = 1.0,       # µm⁻²
    g1 = 30.0,      # µM⁻¹ s⁻¹
    kon2d = 1.0,    # µm² s⁻¹
    kP = 1.0,       # s⁻¹
    nKP = 1.0,
    lambdaX = 1e-4, # µm⁻²
    nC = 3.0,
    X_O1 = 1e-2,    # dimensionless
    O1max = 0.95,
)

# Reasonable fit parameters from an old manual fitting attempt.
const PARAMS_FITTED = (
    fA = 0.0010760218351601953,
    alpha = 149.34436781578844,
    TT = 1.094630372285705,
    g1 = 0.03204443607404042,
    kon2d = 0.08674689927211769,
    kP = 1.6476071717517033,
    nKP = 2.0248129398409676,
    lambdaX = 0.0177085161359101,
    nC = 0.8070635083920701,
    X_O1 = 0.001188267819200772,
    O1max = 0.8571992797935549,
)

# Peptide KDs (µM) from Aleksic et al. 2010.
const PEPTIDE_KD = (var"9V" = 7.97, var"6V" = 22.8, var"3Y" = 47.3, var"4D" = 140.0, var"4A" = 299.0)

"Bound TCR–pMHC density C_T (µm⁻²), the smaller root of the Lauffenburger quadratic."
function get_CT(x, KD, p)
    PT = p.alpha * x
    b = PT + p.TT + p.g1 * KD / p.kon2d
    # This is algebraically equal to (b - sqrt(b² - 4·PT·TT))/2.
    # That form loses precision when PT ≪ b, which happens at low pMHC.
    return p.fA * 2 * PT * p.TT / (b + sqrt(b^2 - 4 * PT * p.TT))
end

"Fully proofread complex density C_N (µm⁻²) after n_KP proofreading steps."
get_CN(CT, KD, p) = (1 / (1 + KD * p.g1 / p.kP))^p.nKP * CT

"Output of the cytoplasmic switch X, in [0, 1]."
get_X(CN, p) = CN^p.nC / (CN^p.nC + p.lambdaX^p.nC)

"Fraction of CD69+ cells O1, in [0, O1max]."
get_O1(X, p) = p.O1max * X / (p.X_O1 + X)

"Fraction of CD69+ cells for pMHC concentration x (µM) and dissociation constant KD (µM)."
function cd69_fraction(x, KD, p)
    CN = get_CN(get_CT(x, KD, p), KD, p)
    return get_O1(get_X(CN, p), p)
end

"Non-degenerate combinations of parameters."
function composite_params(p)
    return (
        xstar = p.TT / p.alpha,                 # µM
        ystar = p.lambdaX / (p.fA * p.TT),      # dimensionless
        KD_Lauf = p.TT * p.kon2d / p.g1,        # µM
        KD_McK = p.kP / p.g1,                   # µM
    )
end

# Discrimination power, as in Pettmann et al. 2021 (eLife 10:e67092) and `get_discrimination_power` in KPXO_model.
# Potency is P15, the pMHC concentration that makes 15% of cells CD69+ (an absolute level, not 15% of each peptide's own top).
# The discrimination power α is the slope of a straight-line fit of log10 P15 against log10 KD, over the peptide panel.

"log10 of the pMHC concentration (µM) at which the response reaches `level`, found by bisection; NaN if it never does."
function log10_potency(KD, p; level = 0.15, lo = -12.0, hi = 4.0)
    cd69_fraction(10.0^hi, KD, p) < level && return NaN
    for _ in 1:60
        mid = (lo + hi) / 2
        cd69_fraction(10.0^mid, KD, p) < level ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end

"""
Fit log10 P15 = intercept + α·log10 KD over the peptides that reach 15%, and return `(; alpha, intercept, log10_p15)`.
Peptides that never reach 15% are left out, as Pettmann et al. leave out potencies that need extrapolation; with fewer than two left, α is NaN.
"""
function discrimination_power(p; KDs = collect(values(PEPTIDE_KD)))
    log10_p15 = [log10_potency(KD, p) for KD in KDs]
    ok = isfinite.(log10_p15)
    sum(ok) < 2 && return (alpha = NaN, intercept = NaN, log10_p15 = log10_p15)
    intercept, alpha = [ones(sum(ok)) log10.(KDs[ok])] \ log10_p15[ok]
    return (; alpha, intercept, log10_p15)
end

"Local slope d log10 P15 / d log10 KD at one KD, by a central difference."
function local_p15_slope(KD, p; h = 1e-3)
    return (log10_potency(KD * 10^h, p) - log10_potency(KD / 10^h, p)) / (2h)
end

# Switching across KD at a fixed pMHC concentration x.
# The response falls from its tight-binding top, cd69_fraction(x, 0, p), toward zero as KD grows.
# KD50 is where it is half that top: the midpoint of the switch in KD.

"log10 of the KD (µM) at which the response at pMHC concentration `x` (µM) is half its tight-binding top, found by bisection; NaN if the top is zero."
function log10_kd50(x, p; lo = -6.0, hi = 8.0)
    half = cd69_fraction(x, 0.0, p) / 2
    half > 0 || return NaN
    for _ in 1:60
        mid = (lo + hi) / 2
        cd69_fraction(x, 10.0^mid, p) > half ? (lo = mid) : (hi = mid)
    end
    return (lo + hi) / 2
end
