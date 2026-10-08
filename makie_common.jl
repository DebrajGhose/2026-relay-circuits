@isdefined(PARULA) || include("common.jl")

import CairoMakie as CM

const MK_THEME = CM.Theme(
    fontsize = 13,
    Axis = (xgridvisible = true, ygridvisible = true,
            xgridcolor = (:black, 0.12), ygridcolor = (:black, 0.12),
            xgridwidth = 0.8, ygridwidth = 0.8,
            rightspinevisible = true, topspinevisible = true,
            spinewidth = 1.0, xlabelsize = 13, ylabelsize = 13,
            xticklabelsize = 11, yticklabelsize = 11),
    Legend = (framevisible = true, framecolor = (:black, 0.5), framewidth = 0.8,
              labelsize = 10, padding = (6, 6, 4, 4), rowgap = 1,
              patchsize = (18, 10)))

function mk_panel_label!(gl, letter; text = "")
    lab = isempty(text) ? letter : letter * "      " * text
    CM.Label(gl, lab; fontsize = 15, halign = :left, tellwidth = false,
             padding = (0, 0, 2, 0))
    return nothing
end

function savemakie_fig(fig, name)
    path = joinpath(FIGDIR, name)

    CM.save(path, fig; px_per_unit = 1.5)
    println("  saved  ", relpath(path, @__DIR__))
    pdfpath = joinpath(FIGDIR, first(splitext(name)) * ".pdf")
    CM.save(pdfpath, fig)
    println("  saved  ", relpath(pdfpath, @__DIR__))
    return path
end

circle_pts(cx, cy, r; n = 240) =
    [CM.Point2f(cx + r*cos(th), cy + r*sin(th)) for th in range(0, 2pi; length = n)]
