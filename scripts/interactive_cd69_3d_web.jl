# Open the CD69 figure in a web view you can rotate with the mouse, using WGLMakie.
# This avoids GLMakie's native window, which lags on some macOS setups (for example on an external monitor).
# Left-drag rotates, scroll zooms, right-drag moves the view.
#
# From a terminal, at the repository root:
#     julia --project=. scripts/interactive_cd69_3d_web.jl
# The figure opens in a browser tab; press Enter in the terminal to quit.
#
# From the VS Code Julia REPL, `include` this file instead.
# The figure appears in the VS Code plot pane, and `fig` and `ax` stay in the REPL.
# After rotating, `ax.azimuth[]` and `ax.elevation[]` give the angles to copy into `VIEWS` in src/plot_cd69_3d.jl.

using WGLMakie

include(joinpath(@__DIR__, "..", "src", "plot_cd69_3d.jl"))

# Loading src/plot_cd69_3d.jl also loads CairoMakie, so choose WGLMakie again here.
WGLMakie.activate!()

fig, ax = cd69_figure(PARAMS_FITTED)
display(fig)

if !isinteractive()
    println("Figure is open in the browser. Press Enter to quit.")
    readline()
end
