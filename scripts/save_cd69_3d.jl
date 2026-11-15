# Write the 3d CD69 figure, with the fitted parameters, as a vector pdf.
#
# From a terminal, at the repository root:
#     julia --project=. scripts/save_cd69_3d.jl
# The figure is written to figures/cd69_3d.pdf, and its path is printed.

include(joinpath(@__DIR__, "..", "src", "plot_cd69_3d.jl"))

println(save_cd69_3d(PARAMS_FITTED))
