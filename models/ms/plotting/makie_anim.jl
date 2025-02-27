using NetCDF, CairoMakie, GLMakie, MathTeXEngine, Statistics

# time-steps
const dt = 0.05
const dto = 100*dt
# space-steps
const dx = 0.1
const dxo = 4*dx
# plotting vectors
const px = collect(dx.*(0:511))
const ox = collect(dxo.*(0:127))

function plot!(ax::Axis, xf::Vector{T}, xa::Vector{T}, xo::Vector{T}, xt::Vector{T}) where {T<:Real}
	plot!(ax, xt, color=:black, label="Truth")
	plot!(ax, xo, color=:black, label="Observations")
	plot!(ax, xf, color = 1, colormap = :tab10, colorrange = (1, 10), label="Forecast")
	plot!(ax, xa, color = 2, colormap = :tab10, colorrange = (1, 10), label="Analysis")
	return nothing
end

function selectRange(t::Vector{T}, tt::Vector{T}) where {T<:Real}
	i = findall(x->any(isapprox.(x, tt; atol = dt*0.5, rtol = 0.0)), t)
	return i
end

function printParameters(p::Vector{T}; parnames = ("\\sigma", "k", "\\tau_{i}", "\\tau_{u}", "\\tau_{o}", "\\tau_{c}", "u_{g}")) where {T<:Real}
	s = ""
	for (i, (par, name)) in enumerate(zip(p, parnames))
		s *= "\\delta " * name * " = " * string(abs(round(par; sigdigits=4))) * ", \\,"
	end
	return s
end

function truthMovie(tt, xt, xo; framerate = 60, name = "truth_anim.mp4")
	GLMakie.activate!()
	# set up animation observables
	n = Observable(1)
	_tt = @lift(round(tt[$n]; digits=2))
	_xt = @lift(xt[1:512, $n])
	_pt = @lift(round.(xt[1025:end, $n]; sigdigits=4))

	_xo = @lift(xo[1:128, $n])

	# set up the figure
	fig = Figure(size = (1024, 512), figure_padding = 8, fontsize = 22, fonts = (; regular = texfont(), bold = texfont()))
	ax = Axis(fig[1:2,1:2], xlabel = L"$x$ [s.u.]", ylabel = L"$u_1(t,x)$")

	lines!(ax, px, _xt, color=:black, label=L"$u_1^{\text{true}}(t, x)$")
	scatter!(ax, ox, _xo, markersize = 6, color=:black, label=L"$u_1^{\text{obs}}(t, x)$")

	xlims!(ax, (minimum(px), maximum(px)))
	ylims!(ax, (-0.1, 1.1))

	axislegend(ax, merge=true, unique=true, orientation = :horizontal)

	tit = Label(fig[0, :], text = L"$t = %$(_tt[]) $ [t.u.]", fontsize = 22, tellwidth=false)

	resize_to_layout!(fig)

	# set animation parameters
	indices = range(0, min(size(xt,2),size(xo,2))-1, step=1)

	# format the output name
	if occursin(".nc", name)
		name = first(split(name, ".nc"))
	elseif occursin(".mp4", name)
		name = first(split(name, ".nc"))
	end
	name = last(split(name, "/"))

	# record animation
	record(fig, "$(name).mp4", indices; framerate = framerate, profile = "high444", compression = 0) do t
		n[] = t+1
		tit.text[] = L"$t = %$(_tt[]) $ [t.u.]"
	end
	return fig, ax
end

function assimMovie(tt, xt, xo, xf, xa; framerate = 60, name = "assim_anim.mp4")
	GLMakie.activate!()
	# set up animation observables
	n = Observable(1)
	_tt = @lift(round(tt[$n]; digits=2))
	_xt = @lift(xt[1:512, $n])
	_pt = @lift(round.(xt[1025:end, $n]; sigdigits=4))

	_xo = @lift(xo[1:128, $n])

	_xf = @lift(xf[1:512, $n])
	_pf = @lift(round.(xf[1025:end, $n]; sigdigits=4))

	_xa = @lift(xa[1:512, $n])
	_pa = @lift(round.(xa[1025:end, $n]; sigdigits=4))

	# set up the figure
	fig = Figure(size = (1024, 512), figure_padding = 8, fontsize = 22, fonts = (; regular = texfont(), bold = texfont()))
	ax = Axis(fig[1:2,1:2], xlabel = L"$x$ [s.u.]", ylabel = L"$u_1(t,x)$")

	lines!(ax, px, _xt, color=:black, label=L"$u_1^{\text{true}}(t, x)$")
	scatter!(ax, ox, _xo, markersize = 6, color=:black, label=L"$u_1^{\text{obs}}(t, x)$")
	lines!(ax, px, _xf, color = 1, colormap = :tab10, colorrange = (1, 10), label=L"$u_1^{\text{b}}(t, x)$")
	lines!(ax, px, _xa, color = 2, colormap = :tab10, colorrange = (1, 10), label=L"$u_1^{\text{a}}(t, x)$")

	xlims!(ax, (minimum(px), maximum(px)))
	ylims!(ax, (-0.1, 1.1))

	axislegend(ax, merge=true, unique=true, orientation = :horizontal)

	tit = Label(fig[0, :], text = L"$t = %$(_tt[]) $ [t.u.]; $ %$(printParameters(_pa[] .- _pt[])) $", fontsize = 22, tellwidth=false)

	resize_to_layout!(fig)

	# set animation parameters
	indices = range(0, min(size(xt,2),size(xo,2),size(xa,2),size(xf,2))-1, step=1)


	# format the output name
	if occursin(".nc", name)
		name = first(split(name, ".nc"))
	elseif occursin(".mp4", name)
		name = first(split(name, ".nc"))
	end
	name = last(split(name, "/"))

	# record animation
	record(fig, "$(name).mp4", indices; framerate = framerate, profile = "high444", compression = 0) do t
		n[] = t+1
		tit.text[] = L"$t = %$(_tt[]) $ [t.u.]; $ %$(printParameters(_pa[] .- _pt[])) $"
		#par.text[] = L"$%$(printParameters(_pa[] .- _pt[]))$"
	end
	return fig, ax
end

function main(ARGS)

	assim_file 	= ARGS[1]
	obs_file 	= ARGS[2]
	truth_file 	= ARGS[3]

	# load the time, forecast, and analysis states from the assim file
	ta = ncread(assim_file, "time")
	xf = ncread(assim_file, "state_for")
	xa = ncread(assim_file, "state_ana")

	# Load the observations from the obs file
	obsmask = map(o->parse(Bool, o), readlines("../obsmask.txt"))
	xo = ncread(obs_file, "obs"); xo = xo[obsmask,:]
	to = ncread(obs_file, "time")

	# Load the truth from the truth file
	xt = ncread(truth_file, "state")
	tt = ncread(truth_file, "time")

	# select the overlapping index ranges
	it = selectRange(tt, ta)
	io = selectRange(to, ta)
	ia = selectRange(ta, to[io])

	# select the overlapping times
	tt = tt[it]
	ta = ta[ia]
	to = to[io]

	# select the overlapping states
	xt = xt[:,it]
	xf = xf[:,ia]
	xa = xa[:,ia]
	xo = xo[:,io]

	# create the truth movie
	truthMovie(tt, xt, xo; name = truth_file);

	# create the assimilation movie
	assimMovie(tt, xt, xo, xf, xa; name = assim_file);

end

main(ARGS)


function plotParametersOverTime(;truth_file = "/home/chris/Development/PDAF/models/ms/state.nc", assim_file = "/home/chris/Development/PDAF/models/ms/FILTER6_FORGET1_FORCING0.01.nc", name = "param_plot.svg", parnames = ("\\sigma", "k", "\\tau_{i}", "\\tau_{u}", "\\tau_{o}", "\\tau_{c}", "u_{g}"))

	# load the time, forecast, and analysis states from the assim file
	ta = ncread(assim_file, "time")
	xa = ncread(assim_file, "state_ana")

	# Load the truth from the truth file
	xt = ncread(truth_file, "state")
	tt = ncread(truth_file, "time")

	# select the overlapping index ranges
	it = selectRange(tt, ta)
	ia = selectRange(ta, tt[it])
#=
	# select the overlapping times
	tt = tt[it]
	ta = ta[ia]

	# select the overlapping states
	xt = xt[:,it]
	xa = xa[:,ia]
=#
	# set up the figure
	CairoMakie.activate!()
	fig = Figure(size = (1024, 512), figure_padding = 8, fontsize = 22, fonts = (; regular = texfont(), bold = texfont()))
	ax = Axis(fig[1:2,1:2], xlabel = L"$t$ [t.u.]", ylabel = L"$\mathbf{p}$", yscale = log10)

	for i in 1:7
		lines!(ax, tt, xt[1024+i,:], color = :black, label=L"$%$(parnames[i])^{\text{true}}$")
	end
	for i in 1:7
		lines!(ax, ta, xa[1024+i,:], color = i, colormap = :tab10, colorrange = (1, 10), linestyle = :dash, linewidth = 2, label=L"$%$(parnames[i])^{\text{a}}$")
		#lines!(ax, tt, xf[1024+i,:], color = i, colormap = :tab10, colorrange = (1, 10), linestyle = :dash, linewidth = 2, label=L"$%$(parnames[i])^{\text{b}}$")
	end

	xlims!(ax, (minimum(tt), maximum(tt)))

	#axislegend(ax, merge=true, unique=true, orientation = :horizontal)
	fig[0, :] = Legend(fig, ax, L"Parameters $ $", framevisible = false, merge=true, unique=true, orientation = :horizontal, ncols = 3)

	resize_to_layout!(fig)

	save(name, fig; pt_per_unit = 1)
	return fig, ax

end
