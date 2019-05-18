#Script for plotting blockade shift results
#26/07/2017

using Plots, JLD, LaTeXStrings
pyplot()
include("functions.jl")
PyPlot.close("all")

atom = "87Rb"
nn = 50
ll = 1
jj = 0.5
mj = 0.5

RRSI, θ, blockadeshiftmeshGHz, C_6val = BlockadeShift(atom,nn,ll,jj,mj)

PyPlot.figure()
heatmap(RRSI*1e6,θ,blockadeshiftmeshGHz*1e3)
plot!(xlabel=L"R\, (μ\mathrm{m})")
plot!(ylabel=L"θ")
plot!(zlabel=L"ΔW\, (\mathrm{MHz})")
plot!(colorbar=false)
plot!(title=string("Calculated Rydberg blockade shift\n for |n,l,j,mⱼ⟩ = |",string(nn),",",string(ll),",",string(jj),",",string(mj),"⟩ state of ",atom))
gui()
# savefig("SaveFiles/BlockadeSurface.pdf")
