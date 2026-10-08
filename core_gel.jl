@isdefined(PARULA) || include("common.jl")
@isdefined(MK_THEME) || include("makie_common.jl")

const CHEMO = (well_radius = 3200.0,
               source_radius = 1000.0,
               cell_diameter = 10.0,
               v0 = 10.0,
               Dr = 0.05,
               sensing_radius = 5.0,
               bias_weight = 0.5,
               sensitivity = 10.0,
               Dc = 500.0,
               kc = 0.026156497379620575,
               Qc = 200_000.0,
               Kd = 20.0,
               dt = 0.5,
               dt_ref = 0.5)

function random_point_in_disc(rng, R)
    rr = R * sqrt(rand(rng)); ph = 2pi*rand(rng)
    return (rr*cos(ph), rr*sin(ph))
end

const GEL = let Rw = CHEMO.well_radius,
                a  = CHEMO.source_radius,
                V  = 80.0 * 1e9,
                lam = sqrt(CHEMO.Dc / CHEMO.kc)
    v_tumour = (2/3) * pi * a^3
    H = (V + v_tumour) / (pi * Rw^2)

    C_surface = 22.368851135413454
    (well_radius = Rw, tumour_radius = a, gel_volume = V, gel_height = H,
     lambda = lam, C_surface = C_surface,
     n_sensor = 5000, n_actuator = 20000,
     bias = CHEMO.bias_weight, hours = 168.0, readout_h = 72.0,
     sensor_sweep = [5, 50, 100, 500, 1000, 2000, 5000, 10000, 20000, 40000],
     reps = 5,

     trace_hours = 168.0, trace_seed_plain = 3, trace_seed_bias = 24)
end

function sensor_ticks(vals, n_actuator)
    labs = [@sprintf("%d\n(%.2f)", vv, vv/(vv + n_actuator)) for vv in vals]
    return (collect(vals), labs)
end

@inline function conc3d(r)
    r <= GEL.tumour_radius && return GEL.C_surface
    return GEL.C_surface * (GEL.tumour_radius/r) * exp(-(r - GEL.tumour_radius)/GEL.lambda)
end

@inline function dfro_bias3d(x, y, z)
    r = sqrt(x*x + y*y + z*z)
    (r <= GEL.tumour_radius || r < 1e-9) && return (0.0, 0.0, 0.0)
    ux, uy, uz = x/r, y/r, z/r
    d = CHEMO.sensing_radius
    fro(c) = c / (c + CHEMO.Kd)
    signal = fro(conc3d(r - d)) - fro(conc3d(r + d))
    mag = CHEMO.sensitivity * signal
    return (-mag*ux, -mag*uy, -mag*uz)
end

@inline function reflect_in_gel(x, y, z, nx, ny, nz)
    Rw = GEL.well_radius; H = GEL.gel_height
    r = sqrt(x*x + y*y)
    if r > Rw
        over = r - Rw; ux, uy = x/r, y/r
        x -= 2*over*ux; y -= 2*over*uy
        dn = nx*ux + ny*uy
        nx -= 2*dn*ux; ny -= 2*dn*uy
    end
    if z < 0.0
        z = -z; nz = -nz
    elseif z > H
        z = 2H - z; nz = -nz
    end
    return (x, y, z, nx, ny, nz)
end

@inline function gel_step(rng, x, y, z, nx, ny, nz, dt, gradient::Bool;
                          bw = GEL.bias)
    s = sqrt(2 * CHEMO.Dr * dt)
    nx += s*randn(rng); ny += s*randn(rng); nz += s*randn(rng)
    m = sqrt(nx*nx + ny*ny + nz*nz)
    nx /= m; ny /= m; nz /= m
    if gradient
        bx, by, bz = dfro_bias3d(x, y, z)
        scale = bw * dt / CHEMO.dt_ref
        nx += scale*bx; ny += scale*by; nz += scale*bz
        m = sqrt(nx*nx + ny*ny + nz*nz)
        m > 1e-12 && (nx /= m; ny /= m; nz /= m)
    end
    x += CHEMO.v0 * dt * nx; y += CHEMO.v0 * dt * ny; z += CHEMO.v0 * dt * nz
    return reflect_in_gel(x, y, z, nx, ny, nz)
end

@inline function closest_approach3(dx0, dy0, dz0, dx1, dy1, dz1)
    ex = dx1 - dx0; ey = dy1 - dy0; ez = dz1 - dz0
    e2 = ex*ex + ey*ey + ez*ez
    e2 < 1e-18 && return dx0*dx0 + dy0*dy0 + dz0*dz0
    tau = -(dx0*ex + dy0*ey + dz0*ez) / e2
    tau = tau < 0 ? 0.0 : (tau > 1 ? 1.0 : tau)
    px = dx0 + tau*ex; py = dy0 + tau*ey; pz = dz0 + tau*ez
    return px*px + py*py + pz*pz
end

@inline function random_polarity(rng)
    uz = 2*rand(rng) - 1; ph = 2pi*rand(rng)
    s = sqrt(max(0.0, 1 - uz*uz))
    return (s*cos(ph), s*sin(ph), uz)
end

@inline function hits_active_sensor3(bi, bj, bk, nb, nbz, starts, items, S0, S,
                                     ax0, ay0, az0, ax, ay, az, sigma2)
    for ii in max(1, bk-1):min(nbz, bk+1)
        for iii in max(1, bj-1):min(nb, bj+1)
            for iiii in max(1, bi-1):min(nb, bi+1)
                bb = ((ii-1)*nb + (iii-1))*nb + iiii
                for iiiii in starts[bb]:(starts[bb+1]-1)
                    ss = items[iiiii]
                    if closest_approach3(S0[1,ss]-ax0, S0[2,ss]-ay0, S0[3,ss]-az0,
                                         S[1,ss]-ax,   S[2,ss]-ay,   S[3,ss]-az) <= sigma2
                        return true
                    end
                end
            end
        end
    end
    return false
end

function run_gel(; gradient::Bool, bw = GEL.bias,
                   n_sensor = GEL.n_sensor, n_actuator = GEL.n_actuator,
                   hours = GEL.hours, dt = CHEMO.dt, seed = 1)
    rng = MersenneTwister(seed)
    Rw = GEL.well_radius; H = GEL.gel_height
    Rt = GEL.tumour_radius; sigma = CHEMO.cell_diameter

    S = zeros(3, n_sensor); SN = zeros(3, n_sensor)
    for ii in 1:n_sensor
        xx, yy = random_point_in_disc(rng, Rw)
        S[1,ii] = xx; S[2,ii] = yy; S[3,ii] = H
        SN[1,ii], SN[2,ii], SN[3,ii] = random_polarity(rng)
    end
    A = zeros(3, n_actuator); AN = zeros(3, n_actuator)
    for ii in 1:n_actuator
        xx, yy = random_point_in_disc(rng, Rw)
        A[1,ii] = xx; A[2,ii] = yy; A[3,ii] = H
        AN[1,ii], AN[2,ii], AN[3,ii] = random_polarity(rng)
    end

    ison = falses(n_actuator)
    ton  = fill(NaN, n_actuator)

    reach = sigma + 2*CHEMO.v0*dt
    binw  = max(100.0, reach)
    half  = Rt + 2*binw
    nb    = ceil(Int, 2*half/binw)
    nbz   = ceil(Int, (Rt + 2*binw)/binw)
    counts = zeros(Int, nb*nb*nbz)
    starts = zeros(Int, nb*nb*nbz + 1)
    items  = zeros(Int, n_sensor)
    active = Int[]; sizehint!(active, n_sensor)

    @inline binof(x, y, z) = begin
        bi = floor(Int, (x + half)/binw) + 1
        bj = floor(Int, (y + half)/binw) + 1
        bk = floor(Int, z/binw) + 1
        (bi < 1 || bi > nb || bj < 1 || bj > nb || bk < 1 || bk > nbz) ? 0 :
            ((bk-1)*nb + (bj-1))*nb + bi
    end

    S0 = similar(S); A0 = similar(A)
    nstep = round(Int, hours*60/dt)
    sigma2 = sigma*sigma
    Rt2 = Rt*Rt
    t = 0.0
    for ii in 1:nstep
        copyto!(S0, S); copyto!(A0, A)
        for iii in 1:n_sensor
            S[1,iii], S[2,iii], S[3,iii], SN[1,iii], SN[2,iii], SN[3,iii] =
                gel_step(rng, S[1,iii], S[2,iii], S[3,iii],
                         SN[1,iii], SN[2,iii], SN[3,iii], dt, gradient; bw = bw)
        end
        for iii in 1:n_actuator
            ison[iii] && continue
            A[1,iii], A[2,iii], A[3,iii], AN[1,iii], AN[2,iii], AN[3,iii] =
                gel_step(rng, A[1,iii], A[2,iii], A[3,iii],
                         AN[1,iii], AN[2,iii], AN[3,iii], dt, gradient; bw = bw)
        end
        t += dt

        empty!(active)
        for iii in 1:n_sensor
            (S[1,iii]^2 + S[2,iii]^2 + S[3,iii]^2 <= Rt2) && push!(active, iii)
        end
        isempty(active) && continue

        fill!(counts, 0)
        for iii in active
            bb = binof(S[1,iii], S[2,iii], S[3,iii]); bb > 0 && (counts[bb] += 1)
        end
        starts[1] = 1
        for iii in 1:length(counts)
            starts[iii+1] = starts[iii] + counts[iii]
        end
        fill!(counts, 0)
        for iii in active
            bb = binof(S[1,iii], S[2,iii], S[3,iii])
            if bb > 0
                items[starts[bb] + counts[bb]] = iii
                counts[bb] += 1
            end
        end

        for iii in 1:n_actuator
            ison[iii] && continue
            bi = floor(Int, (A[1,iii] + half)/binw) + 1
            bj = floor(Int, (A[2,iii] + half)/binw) + 1
            bk = floor(Int, A[3,iii]/binw) + 1
            (bi < 0 || bi > nb+1 || bj < 0 || bj > nb+1 || bk < 0 || bk > nbz+1) && continue
            if hits_active_sensor3(bi, bj, bk, nb, nbz, starts, items, S0, S,
                                   A0[1,iii], A0[2,iii], A0[3,iii],
                                   A[1,iii],  A[2,iii],  A[3,iii], sigma2)
                ison[iii] = true; ton[iii] = t
            end
        end
    end
    return ton, count(isnan, ton)
end

function record_gel_pair(; gradient::Bool, bw = GEL.bias, seed = 1,
                           hours = GEL.trace_hours, keep = 8)
    rng = MersenneTwister(seed)
    Rt = GEL.tumour_radius; sigma = CHEMO.cell_diameter
    r0 = GEL.well_radius - CHEMO.cell_diameter/2
    sx, sy, sz =  r0, 0.0, GEL.gel_height
    ax, ay, az = -r0, 0.0, GEL.gel_height
    snx, sny, snz = random_polarity(rng)
    anx, any_, anz = random_polarity(rng)
    Ps = [sx, sy, sz]; Pa = [ax, ay, az]
    t = 0.0; t_on = NaN
    nstep = round(Int, hours*60/CHEMO.dt)
    for ii in 1:nstep
        sx0, sy0, sz0, ax0, ay0, az0 = sx, sy, sz, ax, ay, az
        sx, sy, sz, snx, sny, snz =
            gel_step(rng, sx, sy, sz, snx, sny, snz, CHEMO.dt, gradient; bw = bw)
        ax, ay, az, anx, any_, anz =
            gel_step(rng, ax, ay, az, anx, any_, anz, CHEMO.dt, gradient; bw = bw)
        t += CHEMO.dt
        if isnan(t_on) && (sx*sx + sy*sy + sz*sz <= Rt*Rt) &&
           closest_approach3(sx0-ax0, sy0-ay0, sz0-az0, sx-ax, sy-ay, sz-az) <= sigma*sigma
            t_on = t
        end
        if ii % keep == 0
            append!(Ps, (sx, sy, sz)); append!(Pa, (ax, ay, az))
        end
    end
    return reshape(Ps, 3, :), reshape(Pa, 3, :), t_on
end

function tumour_occupancy(; gradient::Bool, bw = GEL.bias, hours = GEL.readout_h,
                            n = 500, seed = 1)
    rng = MersenneTwister(seed)
    Rw = GEL.well_radius; H = GEL.gel_height; Rt = GEL.tumour_radius
    P = zeros(3, n); N = zeros(3, n)
    for ii in 1:n
        xx, yy = random_point_in_disc(rng, Rw)
        P[1,ii] = xx; P[2,ii] = yy; P[3,ii] = H
        N[1,ii], N[2,ii], N[3,ii] = random_polarity(rng)
    end
    nstep = round(Int, hours*60/CHEMO.dt)
    for ii in 1:nstep
        for iii in 1:n
            P[1,iii], P[2,iii], P[3,iii], N[1,iii], N[2,iii], N[3,iii] =
                gel_step(rng, P[1,iii], P[2,iii], P[3,iii],
                         N[1,iii], N[2,iii], N[3,iii], CHEMO.dt, gradient; bw = bw)
        end
    end
    inside = count(iii -> P[1,iii]^2 + P[2,iii]^2 + P[3,iii]^2 <= Rt*Rt, 1:n)
    return inside/n, mean(H .- P[3, :])
end

const TRACE_KEEP_3D = 8

function fig_matrigel_3d(; hours = GEL.hours, readout_h = GEL.readout_h,
                           reps = GEL.reps, sweep = GEL.sensor_sweep)
    println("the two-cell circuit in 80 uL of 3D Matrigel")
    Rw = GEL.well_radius; H = GEL.gel_height; Rt = GEL.tumour_radius
    bw = GEL.bias
    ccol = parula_colors(3); scol = ccol[1]; acol = ccol[3]

    @printf("    gel: %.0f uL over a well of radius %.0f um gives a cylinder %.0f um tall\n",
            GEL.gel_volume/1e9, Rw, H)
    @printf("    tumour: hemisphere of radius %.0f um, %.2f%% of the gel by volume\n",
            Rt, 100*(2pi/3*Rt^3)/GEL.gel_volume)
    @printf("    field: C = %.1f nM at the tumour surface, %.1f nM at r = 2000 um, %.2f nM at the far corner\n",
            conc3d(Rt), conc3d(2000.0), conc3d(sqrt(Rw^2 + H^2)))
    dd = dfro_bias3d(2000.0, 0.0, 0.0)
    @printf("    ΔFRO across one cell at r = 2000 um: %.2e ; motility in 3D: v0^2/(6 Dr) = %.0f um^2/min\n",
            sqrt(sum(abs2, dd))/CHEMO.sensitivity, CHEMO.v0^2/(6*CHEMO.Dr))

    Deff = CHEMO.v0^2/(6*CHEMO.Dr)
    @printf("    plain diffusion would take H^2/(2 Deff) = %.0f h to cross the gel\n",
            H^2/(2*Deff)/60)
    for gradient in (false, true)
        fin, dep = tumour_occupancy(; gradient = gradient, hours = readout_h, seed = 5)
        @printf("      %-11s at %.0f h: %.3f of cells inside the tumour, mean depth %.0f um of %.0f\n",
                gradient ? "bias = $bw" : "no gradient", readout_h, fin, dep, H)
    end

    ton_no, nev_no = run_gel(; gradient = false, hours = hours, seed = 11)
    ton_gr, nev_gr = run_gel(; gradient = true, bw = bw, hours = hours, seed = 12)
    fr_no = count(!isnan, ton_no)/GEL.n_actuator
    fr_gr = count(!isnan, ton_gr)/GEL.n_actuator
    @printf("    %d sensors, %d actuators, %.0f h\n", GEL.n_sensor, GEL.n_actuator, hours)
    @printf("      no gradient : %.3f of actuators on  (%d never)\n", fr_no, nev_no)
    @printf("      bias = %g    : %.3f of actuators on  (%d never)\n", bw, fr_gr, nev_gr)

    fmean = Dict{Bool, Vector{Float64}}(); fsem = Dict{Bool, Vector{Float64}}()
    rep_rows = Vector{Any}[]
    for gradient in (false, true)
        mus = Float64[]; sems = Float64[]
        for ns in sweep
            fr = Float64[]
            for ii in 1:reps
                sd = 1000*ii + ns + (gradient ? 7 : 0)
                tn, _ = run_gel(; gradient = gradient, bw = bw, n_sensor = ns,
                                  hours = readout_h, seed = sd)
                non = count(!isnan, tn)
                push!(fr, non/GEL.n_actuator)
                push!(rep_rows, [gradient ? "bias1" : "no_sensing", ns, GEL.n_actuator,
                                 ns/(ns + GEL.n_actuator), ii, sd, non, fr[end]])
            end
            push!(mus, mean(fr)); push!(sems, reps > 1 ? std(fr)/sqrt(reps) : 0.0)
            @printf("      %-12s sensors = %4d : %.3f ± %.3f on at %.0f h\n",
                    gradient ? "bias = $bw" : "no gradient", ns, mus[end], sems[end], readout_h)
        end
        fmean[gradient] = mus; fsem[gradient] = sems
    end

    th = range(0, 2pi; length = 200)
    cxx = Rw .* cos.(th); cyy = Rw .* sin.(th)

    grid_h = range(0, hours; length = 400)
    frac(ts) = [count(<=(gg*60), ts)/GEL.n_actuator for gg in grid_h]
    traces = Dict{Bool, Any}()
    gcol = CM.RGBAf(0.42, 0.42, 0.42, 1.0)
    bcol = PARULA[0.55]

    fig = CM.with_theme(MK_THEME) do
    fig = CM.Figure(size = (1560, 620))

    function gelpanel(col, gradient, letter)
        sd = gradient ? GEL.trace_seed_bias : GEL.trace_seed_plain

        Ps, Pa, t_on = record_gel_pair(; gradient = gradient, bw = bw, seed = sd,
                                         hours = GEL.trace_hours, keep = 1)
        traces[gradient] = (Ps, Pa, t_on, sd)
        gl = fig[1, col] = CM.GridLayout()
        mk_panel_label!(gl[1, 1], letter;
                        text = gradient ? "with chemokine sensing  (bias = $bw)" :
                                          "no chemokine sensing")
        ax = CM.Axis(gl[2, 1]; aspect = CM.DataAspect())
        CM.hidedecorations!(ax); CM.hidespines!(ax)
        CM.limits!(ax, -1.14Rw, 1.14Rw, -1.32Rw, 1.16Rw)

        tcol = gradient ? bcol : CM.RGBAf(0.65, 0.65, 0.65, 1.0)
        CM.lines!(ax, circle_pts(0, 0, Rw); linestyle = :dot, linewidth = 2,
                  color = CM.RGBAf(0.55, 0.55, 0.55, 1.0))
        CM.poly!(ax, circle_pts(0, 0, Rt); color = (tcol, gradient ? 0.30 : 0.16),
                 strokecolor = tcol, strokewidth = 1.4)

        st = 1:TRACE_KEEP_3D:size(Ps, 2)
        CM.lines!(ax, Ps[1, st], Ps[2, st]; color = (scol, 0.85), linewidth = 2.2)
        CM.lines!(ax, Pa[1, st], Pa[2, st]; color = (acol, 0.85), linewidth = 2.2)
        CM.scatter!(ax, [Ps[1, 1], Pa[1, 1]], [Ps[2, 1], Pa[2, 1]];
                    marker = :circle, markersize = 10, color = :white,
                    strokecolor = CM.RGBAf(0.3, 0.3, 0.3, 1.0), strokewidth = 1)
        if !isnan(t_on)
            kk = clamp(round(Int, t_on/CHEMO.dt) + 1, 1, size(Pa, 2))
            CM.scatter!(ax, [Pa[1, kk]], [Pa[2, kk]]; marker = :star5, markersize = 20,
                        color = PARULA[0.95],
                        strokecolor = CM.RGBAf(0.25, 0.25, 0.25, 1.0), strokewidth = 1)
        end
        CM.lines!(ax, [-Rw, -Rw + 1000], [-1.15Rw, -1.15Rw]; linewidth = 3,
                  color = CM.RGBAf(0.3, 0.3, 0.3, 1.0))
        CM.text!(ax, -Rw + 500, -1.23Rw; text = "1 mm", fontsize = 10,
                 align = (:center, :center), color = CM.RGBAf(0.3, 0.3, 0.3, 1.0))
        CM.text!(ax, 0, 1.07Rw; text = "top view, looking down the well", fontsize = 10,
                 align = (:center, :center), color = CM.RGBAf(0.4, 0.4, 0.4, 1.0))

        ents = [CM.PolyElement(color = (tcol, gradient ? 0.30 : 0.16),
                               strokecolor = tcol, strokewidth = 1.4),
                CM.LineElement(color = CM.RGBAf(0.55, 0.55, 0.55, 1.0),
                               linestyle = :dot, linewidth = 2),
                CM.LineElement(color = scol, linewidth = 3),
                CM.LineElement(color = acol, linewidth = 3),
                CM.MarkerElement(marker = :circle, markersize = 10, color = :white,
                                 strokecolor = CM.RGBAf(0.3, 0.3, 0.3, 1.0), strokewidth = 1)]
        labs = [gradient ? "tumour hemisphere, 2 mm across, secreting" :
                           "tumour hemisphere, 2 mm across, no chemokine",
                @sprintf("80 µL gel, 6.4 mm across, %.1f mm deep", H/1000),
                "sensor cell (CD4)", "actuator cell (CD8)", "start, at the gel surface"]
        if !isnan(t_on)
            push!(ents, CM.MarkerElement(marker = :star5, markersize = 14,
                                         color = PARULA[0.95],
                                         strokecolor = CM.RGBAf(0.25, 0.25, 0.25, 1.0),
                                         strokewidth = 1))
            push!(labs, @sprintf("actuator on, %.1f h", t_on/60))
        end
        CM.Legend(gl[2, 1], ents, labs; tellwidth = false, tellheight = false,
                  halign = :right, valign = :bottom, margin = (4, 4, 4, 4))
        CM.rowgap!(gl, 2)
        return ax
    end
    gelpanel(1, false, "A")
    gelpanel(2, true,  "B")

    gr_ = fig[1, 3] = CM.GridLayout()

    ymax = 1.30*max(maximum(frac(ton_gr)), maximum(frac(ton_no)), 0.05)
    mk_panel_label!(gr_[1, 1], "C")
    axC = CM.Axis(gr_[2, 1]; xlabel = "time  (hours)",
                  ylabel = "fraction of actuators switched on",
                  yticks = 0:0.25:1)
    CM.limits!(axC, 0, hours, -0.1, ymax)
    lC1 = CM.lines!(axC, grid_h, frac(ton_no); linewidth = 3, color = gcol)
    lC2 = CM.lines!(axC, grid_h, frac(ton_gr); linewidth = 3, color = bcol)
    CM.lines!(axC, [readout_h, readout_h], [0, ymax]; linewidth = 2, linestyle = :dot,
              color = gcol)
    CM.text!(axC, readout_h + 0.02*hours, 0.45*ymax;
             text = @sprintf("%.0f h,\nread off in D", readout_h), fontsize = 10,
             align = (:left, :center), color = CM.RGBAf(0.35, 0.35, 0.35, 1.0))
    CM.axislegend(axC, [lC1, lC2],
                  [@sprintf("no sensing, %.2f at %.0f h", fr_no, hours),
                   @sprintf("bias = %g, %.2f at %.0f h", bw, fr_gr, hours)];
                  position = :lt, labelsize = 10)

    ymax2 = 1.30*max(maximum(fmean[true]), maximum(fmean[false]), 0.05)
    mk_panel_label!(gr_[3, 1], "D")
    axD = CM.Axis(gr_[4, 1];
                  xlabel = "sensor cells in the gel  (CD4), and their fraction of all cells",
                  ylabel = "fraction of actuators switched on",
                  xscale = log10,
                  xticks = sensor_ticks([5, 50, 500, 5000, 40000],
                                        GEL.n_actuator),
                  yticks = 0:0.25:1)
    CM.limits!(axD, 0.82*minimum(sweep), 1.22*maximum(sweep), -0.1, ymax2)
    hs = CM.Any[]
    for (cc, mu, se) in ((gcol, fmean[false], fsem[false]), (bcol, fmean[true], fsem[true]))
        CM.lines!(axD, sweep, mu; linewidth = 3, color = cc)
        CM.errorbars!(axD, sweep, mu, se; color = cc, whiskerwidth = 8, linewidth = 1.5)
        push!(hs, CM.scatter!(axD, sweep, mu; marker = :circle, markersize = 10,
                              color = cc, strokecolor = :white, strokewidth = 0.8))
    end

    CM.text!(axD, 0.97*maximum(sweep), 0.14*ymax2;
             text = @sprintf("read off at %.0f h, %d actuators, %d gels per point",
                             readout_h, GEL.n_actuator, reps),
             fontsize = 10, align = (:right, :center), color = CM.RGBAf(0.4, 0.4, 0.4, 1.0))
    CM.axislegend(axD, hs, ["no sensing", @sprintf("bias = %g", bw)];
                  position = :lt, labelsize = 10)

    CM.colsize!(fig.layout, 1, CM.Relative(0.31))
    CM.colsize!(fig.layout, 2, CM.Relative(0.31))
    CM.rowgap!(gr_, 2)
    fig
    end
    savemakie_fig(fig, "figure09_matrigel_circuit.png")

    params =
        [["well_radius", Rw, "um", "96-well plate, 6.4 mm across"],
         ["gel_volume", GEL.gel_volume/1e9, "uL", "Matrigel poured over the tumour"],
         ["gel_height", H, "um", "V/(pi R^2); not a free parameter"],
         ["tumour_radius", Rt, "um", "hemisphere on the floor at the centre"],
         ["tumour_volume_fraction", (2pi/3*Rt^3)/GEL.gel_volume, "-", "share of the gel the tumour occupies"],
         ["n_sensor_panelC", GEL.n_sensor, "cells", "figure09_matrigel_circuit: sensor cells (CD4) in the panel C gel"],
         ["n_actuator", GEL.n_actuator, "cells", "figure09_matrigel_circuit: actuator cells (CD8), fixed throughout"],
         ["sensor_sweep", join(sweep, ";"), "cells", "figure09_matrigel_circuit panel D x axis, drawn on a log scale; the assay's own 0.5k-40k"],
         ["bias_weight", bw, "-", "steering strength with chemokine sensing; 0 without"],
         ["contact_distance", CHEMO.cell_diameter, "um", "an actuator switches on within this of an active sensor"],
         ["v0", CHEMO.v0, "um/min", "crawling speed, fixed"],
         ["Dr", CHEMO.Dr, "rad^2/min", "rotational diffusion of the polarity vector"],
         ["D_eff_3d", CHEMO.v0^2/(6*CHEMO.Dr), "um^2/min", "v0^2/(6 Dr)"],
         ["sensing_radius", CHEMO.sensing_radius, "um", "front to back sensing half distance"],
         ["sensitivity", CHEMO.sensitivity, "-", "receptor sensitivity"],
         ["Dc", CHEMO.Dc, "um^2/min", "chemokine diffusion coefficient"],
         ["kc", CHEMO.kc, "1/min", "CXCL12 removal rate, measured: ln(2)/26.5 min half-life"],
         ["lambda", GEL.lambda, "um", "range of the field, sqrt(Dc/kc)"],
         ["C_surface", GEL.C_surface, "nM", "field at the tumour surface, from the tumour CXCL12 ELISA"],
         ["Kd", CHEMO.Kd, "nM", "receptor dissociation constant"],
         ["dt", CHEMO.dt, "min", "time step"],
         ["hours_panelC", hours, "h", "how long each gel is watched"],
         ["readout_h_panelD", readout_h, "h", "assay time panel D is read off at"],
         ["reps_panelD", reps, "gels", "independent gels per point in panel D"],
         ["trace_hours", GEL.trace_hours, "h", "window drawn in panels A and B"],
         ["trace_seed_plain", GEL.trace_seed_plain, "-", "seed of the panel A run"],
         ["trace_seed_bias", GEL.trace_seed_bias, "-", "seed of the panel B run"],
         ["start", "uniform on the free surface of the gel", "-", "panels C and D; panels A and B start antipodal on the rim"],
         ["field", "C(r) = C_surface (a/r) exp(-(r-a)/lambda)", "-", "3D screened kernel; flat inside the tumour"],
         ["geometry", "3D cylinder of gel over a hemisphere", "-", "gel cylinder over a tumour hemisphere"]]

    data_rows = Vector{Any}[]
    sfrac(ns) = ns/(ns + GEL.n_actuator)
    figname = "figure09_matrigel_circuit"

    for (gradient, panel) in ((false, "A"), (true, "B"))
        Ps, Pa, t_on, sd = traces[gradient]
        cond = gradient ? "bias1" : "no_sensing"
        for (role, P) in (("sensor_CD4", Ps), ("actuator_CD8", Pa))
            for kk in 1:size(P, 2)
                push!(data_rows, figdata_row(; figure = figname,
                      record = "trajectory", panel = panel,
                      condition = cond, cell_role = role, n_sensor = 1, n_actuator = 1,
                      sensor_fraction = 0.5, seed = sd, index = kk,
                      time_h = (kk-1)*CHEMO.dt/60,
                      x_um = P[1,kk], y_um = P[2,kk], z_um = P[3,kk]))
            end
        end
        isnan(t_on) || push!(data_rows, figdata_row(; figure = figname,
              record = "trace_switch",
              panel = panel, condition = cond, cell_role = "actuator_CD8",
              n_sensor = 1, n_actuator = 1, sensor_fraction = 0.5, seed = sd,
              time_h = t_on/60, n_switched_on = 1, frac_on = 1))
    end

    for (gradient, ton, sd) in ((false, ton_no, 11), (true, ton_gr, 12))
        cond = gradient ? "bias1" : "no_sensing"
        for ii in eachindex(ton)
            push!(data_rows, figdata_row(; figure = figname,
                  record = "switch_time", panel = "C",
                  condition = cond, cell_role = "actuator_CD8",
                  n_sensor = GEL.n_sensor, n_actuator = GEL.n_actuator,
                  sensor_fraction = sfrac(GEL.n_sensor), seed = sd, index = ii,
                  time_h = isnan(ton[ii]) ? "" : ton[ii]/60,
                  n_switched_on = isnan(ton[ii]) ? 0 : 1))
        end
    end

    for rr in rep_rows
        cond, ns, nact, sf, rep, sd, non, fr = rr
        push!(data_rows, figdata_row(; figure = figname,
              record = "replicate", panel = "D",
              condition = cond, cell_role = "actuator_CD8", n_sensor = ns,
              n_actuator = nact, sensor_fraction = sf, replicate = rep, seed = sd,
              time_h = readout_h, n_switched_on = non, frac_on = fr))
    end

    return fig, data_rows, params
end
