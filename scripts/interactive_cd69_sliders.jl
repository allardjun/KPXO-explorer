# Show the 3d CD69 figure with a slider for each model parameter, using WGLMakie.
# Moving a slider redraws the five lines; left-drag on the plot still rotates it.
#
# From a terminal, at the repository root:
#     julia --project=. scripts/interactive_cd69_sliders.jl
# The figure opens in a browser tab; press Enter in the terminal to quit.
#
# From the VS Code Julia REPL, `include` this file instead.
# The figure appears in the VS Code plot pane, and `fig`, `ax` and `params` stay in the REPL.
# `params[]` is the parameter set the sliders show now.
#
# The sliders only work while Julia is running, since Julia computes the lines.

using WGLMakie

include(joinpath(@__DIR__, "..", "src", "cd69_sliders.jl"))

# Loading src/plot_cd69_3d.jl also loads CairoMakie, so choose WGLMakie again here.
WGLMakie.activate!()

fig, ax, params = cd69_slider_figure()
display(fig)

if !isinteractive()
    println("Figure is open in the browser. Press Enter to quit.")
    readline()
end
