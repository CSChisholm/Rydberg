#Script to plot Rydberg radial wave functions
#23/07/2017

using Plots, JLD, LaTeXStrings
pyplot()
include("functions.jl")
PyPlot.close("all")

#Input information
atom = "87Rb"
nn = 50
ll = 0
jj = 0.5

#Calculate wave function
normY_sol, rr = numerovfunc(atom,nn,ll,jj)

#Rescale for plotting
plotscale = sqrt(rr)
probamp = (normY_sol.*plotscale).^2

alpha_c = getalpha(atom)

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
plot!(ylabel=L"|rR(r)|^2 \,(a_0^{-1})")
plot!(title=string(atom," radial probability density |n,l,j⟩ = |",string(nn),",",string(ll),",",string(jj),"⟩"))
plot!(leg=false)
gui()
