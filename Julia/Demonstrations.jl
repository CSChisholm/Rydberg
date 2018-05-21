#Script to make a series of demonstration plots
#27/07/2017

using Plots, JLD, LaTeXStrings
pyplot()
include("functions.jl")
PyPlot.close("all")

#Input information
atom = "87Rb"
nn = 50
ll = 0
jj = 0.5
mj = 0.5

alpha_c = getalpha(atom)

#Calculate wave function
normY_sol, rr = numerovfunc(atom,nn,ll,jj)

#Rescale for plotting
plotscale = sqrt(rr)
probamp = (normY_sol.*plotscale).^2

PyPlot.figure()
if nn>20
  plot(plotscale[plotscale .> sqrt(alpha_c^(1/3))],normY_sol[plotscale .> sqrt(alpha_c^(1/3))],linewidth=2)
else
  plot(plotscale,normY_sol,linewidth=2)
end
plot!(xlabel=L"(r/a_0)^{1/2}")
plot!(ylabel=L"r^{1/2}R(r) \,(a_0^{-1})")
plot!(title=string(atom," radial wavefunction |n,l,j⟩ = |",string(nn),",",string(ll),",",string(jj),"⟩"))
plot!(leg=false)
gui()

PyPlot.figure()
if nn>20
  plot(rr[rr .> alpha_c^(1/3)],probamp[rr .> alpha_c^(1/3)],linewidth=2)
else
  plot(rr,probamp,linewidth=2)
end
plot!(xlabel=L"r/a_0")
plot!(ylabel=L"|r^{1/2}R(r)|^2 \,(a_0^{-1})")
plot!(title=string(atom," radial probability density |n,l,j⟩ = |",string(nn),",",string(ll),",",string(jj),"⟩"))
plot!(leg=false)
gui()

#Compare groundstate to self-consistent field method
scfdict = load("Data/SCFcalc.jld")
plotscale2 = scfdict["plotscale2"]
sqrtrR = scfdict["sqrtrR"]
normY_solg, rrg = numerovfunc(atom,5,0,0.5)

PyPlot.figure()
plot(sqrt(rrg),-normY_solg,linewidth=2,label="Numerov")
scatter!(plotscale2,sqrtrR,label="Callaway")
plot!(xlabel=L"r/a_0")
plot!(ylabel=L"|r^{1/2}R(r)|^2 \,(a_0^{-1})")
plot!(title="Rubidium Groundstate Radial Wavefunction")
gui()

#Lifetimes
ltdict = load("Data/Lifetimes.jld")
sexpt = ltdict["sexpt"]
sexpterr = ltdict["sexpterr"]
pexpt = ltdict["pexpt"]
pexpterr = ltdict["pexpterr"]
dexpt = ltdict["dexpt"]
dexpterr = ltdict["dexpterr"]

#S1_2
nns = collect(28:45)
calcs = zeros(length(nns))
for kk = 1:length(nns)
  calcs[kk] = Radiative_lifetimes(atom,nns[kk],0,0.5)
end
calcs = calcs*1e6

#P3_2
nnp = collect(34:44)
calcp = zeros(length(nnp))
for kk = 1:length(nnp)
  calcp[kk] = Radiative_lifetimes(atom,nnp[kk],1,1.5)
end
calcp = calcp*1e6

#D5_2
nnd = collect(29:44)
calcd = zeros(length(nnd))
for kk = 1:length(nnd)
  calcd[kk] = Radiative_lifetimes(atom,nnd[kk],2,2.5)
end
calcd = calcd*1e6

PyPlot.figure()
scatter(nns,calcs,label="Calculated")
scatter!(nns,sexpt[:],yerr=repmat(sexpterr[:],1,2),label="Experimental")
plot!(xlabel=L"n")
plot!(ylabel=L"τ\, (μ\mathrm{s})")
plot!(title="Plot of calculated and experimental lifetimes for \n |n,l,j⟩=|n,0,0.5⟩ states")
gui()

PyPlot.figure()
scatter(nnp,calcp,label="Calculated")
scatter!(nnp,pexpt[:],yerr=repmat(pexpterr[:],1,2),label="Experimental")
plot!(xlabel=L"n")
plot!(ylabel=L"τ\, (μ\mathrm{s})")
plot!(title="Plot of calculated and experimental lifetimes for \n |n,l,j⟩=|n,1,1.5⟩ states")
gui()

PyPlot.figure()
scatter(nnd,calcd,label="Calculated")
scatter!(nnd,dexpt[:],yerr=repmat(dexpterr[:],1,2),label="Experimental")
plot!(xlabel=L"n")
plot!(ylabel=L"τ\, (μ\mathrm{s})")
plot!(title="Plot of calculated and experimental lifetimes for \n |n,l,j⟩=|n,2,2.5⟩ states")
gui()

#Blockade shift
RRSI, θ, blockadeshiftmeshGHz, C_6val = BlockadeShift(atom,nn,ll,jj,mj)

PyPlot.figure()
heatmap(RRSI*1e6,θ,blockadeshiftmeshGHz*1e3)
plot!(xlabel=L"R\, (μ\mathrm{m})")
plot!(ylabel=L"θ")
plot!(zlabel=L"ΔW\, (\mathrm{MHz})")
plot!(title=string("Calculated Rydberg blockade shift for \n |n,l,j,mⱼ⟩ = |",string(nn),",",string(ll),",",string(jj),",",string(mj),"⟩ state of ",atom," [MHz]"))
gui()

#C₆
c6dict = load("Data/C6data.jld")
nnc6 = c6dict["nnc6"]
#Convert from GHz⋅m⁻⁶ to GHz⋅μm⁻⁶
c6S1_2 = c6dict["c6S1_2"]*1e36
c6D3_2 = c6dict["c6D3_2"]*1e36
c6D5_2 = c6dict["c6D5_2"]*1e36

PyPlot.figure()
plot(nnc6,abs(c6S1_2),linewidth=2,marker=(:circle))
plot!(yaxis=(:log))
plot!(xlabel=L"n")
plot!(ylabel=L"|C₆|\, (\mathrm{GHz}⋅μ\mathrm{m}⁻⁶)")
plot!(title="Calculated C₆ for |n,l,j⟩ = |n,0,0.5⟩")
plot!(leg=false)
gui()

PyPlot.figure()
plot(nnc6,abs(c6D3_2),linewidth=2,marker=(:circle))
plot!(yaxis=(:log))
plot!(xlabel=L"n")
plot!(ylabel=L"|C₆|\, (\mathrm{GHz}⋅μ\mathrm{m}⁻⁶)")
plot!(title="Calculated C₆ for |n,l,j⟩ = |n,2,1.5⟩")
plot!(leg=false)
gui()

PyPlot.figure()
plot(nnc6,abs(c6D5_2),linewidth=2,marker=(:circle))
plot!(yaxis=(:log))
plot!(xlabel=L"n")
plot!(ylabel=L"|C₆|\, (\mathrm{GHz}⋅μ\mathrm{m}⁻⁶)")
plot!(title="Calculated C₆ for |n,l,j⟩ = |n,2,2.5⟩")
plot!(leg=false)
gui()
