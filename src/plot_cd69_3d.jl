# Forward model: fraction of CD69+ cells versus pMHC concentration and peptide KD.
# Draws one dose-response line per peptide KD, placed in 3d axes, and faint lines across KD at fixed pMHC concentrations.
#
# This file only defines the plotting functions; scripts/save_cd69_3d.jl writes the figure to figures/cd69_3d.pdf.

using CairoMakie
using Printf

include(joinpath(@__DIR__, "model.jl"))

# Generated figures go in figures/ at the repository root, which git ignores.
const FIGDIR = normpath(joinpath(@__DIR__, "..", "figures"))

# Tick labels as powers of ten, since the axes hold log10 values.
pow10_labels(values) = [rich("10", superscript(string(round(Int, v)))) for v in values]

# Camera angles for the 3d plot, and the axis each one looks straight along, whose ticks are then hidden.
# Makie's Axis3 has no perspective by default, so the side views show true slopes.
# The small tilts keep the ticks on the near edges; an exact 0 puts them on the far edge.
const VIEWS = (
    Default = (azimuth = -0.35π, elevation = 0.15π, along = nothing),
    Dose = (azimuth = -0.5π, elevation = 0.002, along = :y),       # response against pMHC, looking along KD
    KD = (azimuth = 0.0, elevation = 0.002, along = :x),           # response against KD, looking along pMHC
    Top = (azimuth = -0.5π, elevation = 0.5π, along = :z),         # the floor: P15 against KD
)

"Turn the axis to one of `VIEWS`, and hide the ticks and label of the axis it looks along."
function set_view!(ax, view::Symbol)
    v = VIEWS[view]
    ax.azimuth[] = v.azimuth
    ax.elevation[] = v.elevation
    show_axis_decorations!(ax; hide = v.along)
    return ax
end

"""
Show the ticks and labels of all three axes, except the one named by `hide`.
In a side or top view the box fills its space, so the right-hand axis label needs extra room (left, right, bottom, top).
"""
function show_axis_decorations!(ax; hide = nothing)
    for a in (:x, :y, :z)
        visible = a != hide
        for attr in (:labelvisible, :ticklabelsvisible, :ticksvisible)
            getproperty(ax, Symbol(a, attr))[] = visible
        end
    end
    ax.protrusions[] = hide === nothing ? 30 : (30, 80, 30, 30)
end

"One color per peptide, tightest first. Stops short of the yellow end of viridis, which is hard to see on white."
peptide_colors(n) = reverse(get(cgrad(:viridis), range(0, 0.85, length = n)))

# The default pMHC range runs past the tested range (1.0×10⁻⁸ to 0.1 µM, in fivefold steps), so that every curve reaches its plateau.
# The weakest peptide, 4A, reaches 99% of its plateau near 10^1.1 µM.
"""
Draw the 3d plot in a new axis at figure position `pos`, and return the axis.
`p` is a parameter set, or an `Observable` of one; with an `Observable`, the lines follow its changes.
`fixed_log10x` gives the log10 pMHC concentrations of the faint lines across KD; pass an empty range to leave them out.
`potency_floor` draws the P15 marks and the straight-line fit on the floor; they are worth leaving out in a side view, where the floor collapses onto an axis.
`view` is one of the camera angles in `VIEWS`.
"""
function cd69_axis!(pos, p; KDs = collect(values(PEPTIDE_KD)), labels = string.(keys(PEPTIDE_KD)),
                    log10x_range = (-8.0, 2.0), fixed_log10x = log10x_range[1]:1:log10x_range[2],
                    potency_floor = true, view = :Default)
    p = p isa Observable ? p : Observable(p)
    log10x = range(log10x_range..., length = 400)
    log10KD = log10.(KDs)
    log10KD_range = (floor(minimum(log10KD)), ceil(maximum(log10KD)))

    ax = Axis3(pos;
        xlabel = "pMHC concentration (µM)",
        ylabel = "Peptide KD (µM)",
        zlabel = "Fraction CD69+",
        # From the default angle, pMHC runs along the front edge and KD recedes to the back right.
        azimuth = VIEWS.Default.azimuth, elevation = VIEWS.Default.elevation,
        xticks = (log10x_range[1]:2:log10x_range[2], pow10_labels(log10x_range[1]:2:log10x_range[2])),
        yticks = (log10KD_range[1]:log10KD_range[2], pow10_labels(log10KD_range[1]:log10KD_range[2])),
        # The zero tick is left off, since it collides with the top KD tick; the floor marks zero.
        zticks = 0.25:0.25:1,
        limits = (log10x_range, log10KD_range, (0, 1)),
        xlabeloffset = 45, ylabeloffset = 45, zlabeloffset = 50,
    )

    colors = peptide_colors(length(KDs))

    # Faint lines at fixed pMHC concentrations, one per decade, running across a continuous range of KD.
    # Seen with the KD axis horizontal, these show how sharply the response switches with KD.
    KD_fine = 10 .^ range(log10KD_range..., length = 200)
    for lx in fixed_log10x
        fixed_x_line = lift(p -> [Point3f(lx, log10(KD), cd69_fraction(10.0^lx, KD, p)) for KD in KD_fine], p)
        lines!(ax, fixed_x_line; color = (:gray50, 0.6), linewidth = 1)
    end

    # A dot on each of those lines at its KD50, where the response is half its tight-binding top.
    # Dots that fall outside the KD axis are left out (NaN).
    kd50_points = lift(p) do p
        map(fixed_log10x) do lx
            y = log10_kd50(10.0^lx, p)
            inside = log10KD_range[1] <= y <= log10KD_range[2]
            inside ? Point3f(lx, y, cd69_fraction(10.0^lx, 0.0, p) / 2) : Point3f(NaN, NaN, NaN)
        end
    end
    scatter!(ax, kd50_points; color = :gray25, markersize = 6)

    dp = lift(discrimination_power, p)

    if potency_floor
        # Potency P15 on the floor, over a continuous range of KD; gaps where the response never reaches 15%.
        p15_path = lift(p) do p
            [Point3f(log10_potency(KD, p), log10(KD), 0) for KD in KD_fine]
        end
        lines!(ax, p15_path; color = :gray40, linestyle = :dash, linewidth = 1)

        # Straight-line fit through the peptides' P15 values, over the panel's KD range; its slope is the discrimination power.
        fit_ends = log10.(extrema(KDs))
        fit_line = lift(dp -> [Point3f(dp.intercept + dp.alpha * y, y, 0) for y in fit_ends], dp)
        lines!(ax, fit_line; color = :firebrick, linewidth = 1.5)
        alpha_position = lift(fit_line -> fit_line[end], fit_line)
        text!(ax, alpha_position; text = lift(dp -> @sprintf("slope %.2f", dp.alpha), dp), color = :firebrick,
            align = (:left, :center), offset = (8, 0), fontsize = 12)
    end

    for (i, KD) in enumerate(KDs)
        y = fill(log10(KD), length(log10x))
        z = lift(p -> cd69_fraction.(10 .^ log10x, KD, Ref(p)), p)
        # Faint footprint on the floor helps the eye place each line in depth.
        lines!(ax, log10x, y, zeros(length(log10x)); color = (colors[i], 0.35), linewidth = 1)
        lines!(ax, log10x, y, z; color = colors[i], linewidth = 2.5)
        # Labels go at the flat left ends, where the lines are well apart in depth.
        label_position = lift(z -> Point3f(log10x[1], log10(KD), z[1]), z)
        text!(ax, label_position; text = labels[i], color = colors[i],
            align = (:right, :center), offset = (-6, 0), fontsize = 12)
        # This peptide's P15 on the floor.
        potency_floor && scatter!(ax, lift(dp -> Point3f(dp.log10_p15[i], log10(KD), 0), dp); color = colors[i], markersize = 7)
    end

    view == :Default || set_view!(ax, view)
    return ax
end

"Build the 3d figure and return it with its axis, without saving, so any Makie backend can show it."
function cd69_figure(p; kwargs...)
    fig = Figure(size = (620, 520), fontsize = 13)
    ax = cd69_axis!(fig[1, 1], p; kwargs...)
    return fig, ax
end

"Save the figure as a vector pdf in the figures directory and return its path."
function save_cd69_3d(p; filename = "cd69_3d.pdf", kwargs...)
    fig, _ = cd69_figure(p; kwargs...)
    mkpath(FIGDIR)
    path = joinpath(FIGDIR, filename)
    CairoMakie.save(path, fig; backend = CairoMakie)
    return path
end
