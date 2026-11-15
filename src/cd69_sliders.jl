# The 3d CD69 figure with a slider for each model parameter.
# This file only defines `cd69_slider_figure`; scripts/interactive_cd69_sliders.jl shows it.

using Printf

include(joinpath(@__DIR__, "plot_cd69_3d.jl"))

# One slider per parameter: name, label, scale, and bounds.
# A `:log` slider moves in log10 space; the hill exponents and the two fractions move in plain space.
# The bounds cover the nominal and fitted sets.
const SLIDERS = [
    (name = :fA,      label = "fA",             scale = :lin, bounds = (0.0, 1.0)),
    (name = :alpha,   label = "α (µm⁻²/µM)",    scale = :log, bounds = (1e-1, 1e3)),
    (name = :TT,      label = "TT (µm⁻²)",      scale = :log, bounds = (1e-1, 1e3)),
    (name = :g1,      label = "g̃1 (µM⁻¹s⁻¹)",  scale = :log, bounds = (1e-3, 1e2)),
    (name = :kon2d,   label = "kon2d (µm²s⁻¹)", scale = :log, bounds = (1e-3, 1e1)),
    (name = :kP,      label = "kP (s⁻¹)",       scale = :log, bounds = (1e-2, 1e2)),
    (name = :nKP,     label = "nKP",            scale = :lin, bounds = (0.1, 6.0)),
    (name = :lambdaX, label = "λX (µm⁻²)",      scale = :log, bounds = (1e-6, 1e1)),
    (name = :nC,      label = "nC",             scale = :lin, bounds = (0.1, 5.0)),
    (name = :X_O1,    label = "X_O1",           scale = :log, bounds = (1e-5, 1e1)),
    (name = :O1max,   label = "O1max",          scale = :lin, bounds = (0.0, 1.0)),
]

# Values a slider can take, in its own space (log10 for `:log`).
# The exact values of the given parameter sets are added, so the buttons can return to them exactly.
function slider_range(s, psets; steps = 500)
    lo, hi = to_slider.(Ref(s), s.bounds)
    return sort(unique([range(lo, hi, length = steps + 1); [to_slider(s, getfield(p, s.name)) for p in psets]]))
end

to_param(s, v) = s.scale == :log ? 10.0^v : v
to_slider(s, value) = s.scale == :log ? log10(value) : value

"""
Text for the four combinations of parameters that the output depends on, and for the discrimination power.
The local slopes are taken at the tightest and weakest peptides.
"""
function composite_text(p)
    c = composite_params(p)
    dp = discrimination_power(p)
    KD_low, KD_high = extrema(values(PEPTIDE_KD))
    return join([
        @sprintf("x* = TT/α = %.3g µM", c.xstar),
        @sprintf("y* = λX/(fA·TT) = %.3g", c.ystar),
        @sprintf("KD Lauf = TT·kon2d/g̃1 = %.3g µM", c.KD_Lauf),
        @sprintf("KD McK = kP/g̃1 = %.3g µM", c.KD_McK),
        @sprintf("Discrimination power (P15 fit) = %.2f", dp.alpha),
        @sprintf("Local P15 slope: %.2f at %g µM, %.2f at %g µM",
            local_p15_slope(KD_low, p), KD_low, local_p15_slope(KD_high, p), KD_high),
    ], "\n")
end

"Build the figure with sliders; return it with the 3d axis and an `Observable` of the current parameter set."
function cd69_slider_figure(p0 = PARAMS_FITTED; presets = (Fitted = PARAMS_FITTED, Nominal = PARAMS_NOMINAL))
    fig = Figure(size = (1050, 560), fontsize = 13)
    panel = GridLayout(fig[1, 2]; tellheight = false, valign = :center)

    sg = SliderGrid(panel[1, 1],
        [(label = s.label,
          range = slider_range(s, values(presets)),
          format = v -> @sprintf("%.3g", to_param(s, v)),
          startvalue = to_slider(s, getfield(p0, s.name)))
         for s in SLIDERS]...;
        width = 360)

    names = Tuple(s.name for s in SLIDERS)
    params = lift((vs...) -> NamedTuple{names}(Tuple(map(to_param, SLIDERS, vs))),
                  (sl.value for sl in sg.sliders)...)

    ax = cd69_axis!(fig[1, 1], params)

    Label(panel[2, 1], lift(composite_text, params); justification = :left, halign = :left)

    buttons = GridLayout(panel[3, 1]; halign = :left)
    for (i, (preset_name, preset)) in enumerate(pairs(presets))
        b = Button(buttons[1, i]; label = string(preset_name))
        on(b.clicks) do _
            for (s, sl) in zip(SLIDERS, sg.sliders)
                set_close_to!(sl, to_slider(s, getfield(preset, s.name)))
            end
        end
    end
    print_button = Button(buttons[1, length(presets) + 1]; label = "Print to REPL")
    on(_ -> println(params[]), print_button.clicks)

    # Buttons that turn the plot to a side or top view; "KD" shows the response against KD.
    views = GridLayout(panel[4, 1]; halign = :left)
    Label(views[1, 1], "View:")
    for (i, view) in enumerate(keys(VIEWS))
        b = Button(views[1, i + 1]; label = string(view))
        on(_ -> set_view!(ax, view), b.clicks)
    end
    # Rotating away from a side view by hand brings back the hidden ticks.
    onany(ax.azimuth, ax.elevation) do az, el
        current = findfirst(v -> isapprox(v.azimuth, az; atol = 1e-6) && isapprox(v.elevation, el; atol = 1e-6), VIEWS)
        current === nothing && show_axis_decorations!(ax)
    end

    colsize!(fig.layout, 1, Relative(0.6))
    colgap!(fig.layout, 1, 30)
    return fig, ax, params
end
