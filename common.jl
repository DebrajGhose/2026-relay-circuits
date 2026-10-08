import Pkg
Pkg.activate(@__DIR__)

using Statistics, Random, LinearAlgebra, Printf
import CairoMakie as CM

const PARULA = CM.cgrad([
    CM.RGB(0.2081,0.1663,0.5292), CM.RGB(0.0589,0.3678,0.8843), CM.RGB(0.0779,0.5040,0.8384),
    CM.RGB(0.0200,0.5965,0.8154), CM.RGB(0.0723,0.6859,0.6986), CM.RGB(0.4668,0.7271,0.4707),
    CM.RGB(0.7900,0.7250,0.3300), CM.RGB(0.9856,0.7572,0.2502), CM.RGB(0.9763,0.9831,0.0538)])

parula_colors(k) = k == 1 ? [PARULA[0.18]] : [PARULA[tt] for tt in range(0.0, 0.86; length = k)]

const FIGDIR = normpath(joinpath(@__DIR__, "figures"))
isdir(FIGDIR) || mkpath(FIGDIR)

csv_cell(v) = begin
    s = v isa Real ? @sprintf("%.6g", v) : string(v)
    (occursin(',', s) || occursin('"', s)) ? '"' * replace(s, '"' => "\"\"") * '"' : s
end

csv_row(vals) = join((csv_cell(vv) for vv in vals), ",")

const FIGDATA_HEADER = ["figure", "record", "panel", "condition", "cell_role",
                        "n_sensor", "n_actuator", "sensor_fraction_of_all_cells",
                        "replicate", "seed", "index", "time_h",
                        "x_um", "y_um", "z_um", "n_switched_on", "frac_on"]

function figdata_row(; figure, record, panel, condition = "", cell_role = "",
                       n_sensor = "", n_actuator = "", sensor_fraction = "",
                       replicate = "", seed = "", index = "", time_h = "",
                       x_um = "", y_um = "", z_um = "", n_switched_on = "", frac_on = "")
    return Any[figure, record, panel, condition, cell_role, n_sensor, n_actuator,
               sensor_fraction, replicate, seed, index, time_h,
               x_um, y_um, z_um, n_switched_on, frac_on]
end

function savetable_fig(name, header, rows; meta = nothing)
    path = joinpath(FIGDIR, name)
    open(path, "w") do io
        if meta !== nothing
            println(io, "# ", name)
            println(io, "# every parameter this table was produced under follows.")
            println(io, "# parameter,value,units,note")
            for rr in meta
                println(io, "# ", csv_row(rr))
            end
            println(io, "#")
        end
        println(io, csv_row(header))
        for rr in rows
            length(rr) == length(header) ||
                error("savetable_fig: row of $(length(rr)) fields for $(length(header)) columns")
            println(io, csv_row(rr))
        end
    end
    println("  saved  ", relpath(path, @__DIR__))
    return path
end
