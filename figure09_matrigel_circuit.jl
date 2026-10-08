@isdefined(fig_matrigel_3d) || include("core_gel.jl")

function run_figure09()
    fig, rows, params = fig_matrigel_3d()
    savetable_fig("figure09_matrigel_circuit_data.csv", FIGDATA_HEADER, rows; meta = params)
    return fig, rows, params
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_figure09()
end
