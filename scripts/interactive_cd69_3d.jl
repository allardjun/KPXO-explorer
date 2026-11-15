# Open the CD69 figure in a window you can rotate with the mouse.
# Left-drag rotates, scroll zooms, right-drag moves the view.
#
# From a terminal, at the repository root:
#     julia --project=. scripts/interactive_cd69_3d.jl
# The script waits until you close the window.
#
# From the VS Code Julia REPL, `include` this file instead.
# It returns right away and leaves `fig` and `ax` in the REPL.
# After rotating, `ax.azimuth[]` and `ax.elevation[]` give the angles to copy into `VIEWS` in src/plot_cd69_3d.jl.

using GLMakie

include(joinpath(@__DIR__, "..", "src", "plot_cd69_3d.jl"))

# Loading src/plot_cd69_3d.jl also loads CairoMakie, so choose GLMakie again here.
# `inline = false` opens a separate window, instead of a still image in the VS Code plot pane.
# GLMakie's default of 30 frames per second makes dragging feel sluggish; a frame of this figure takes about 8 ms, so 60 is easy.
GLMakie.activate!(inline = false, framerate = 60)

fig, ax = cd69_figure(PARAMS_FITTED)
screen = display(fig)

if !isinteractive()
    wait(screen)
end
