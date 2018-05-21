#Script for calculating blockade shift coefficients
#26/07/2017

using JLD
include("functions.jl")

nnc6 = collect(35:82)
atom = "87Rb"

c6S1_2 = zeros(length(nnc6))
c6D3_2 = zeros(length(nnc6))
c6D5_2 = zeros(length(nnc6))

for kk = 1:length(nnc6)
  junk1, junk2, junk3, c6S1_2[kk] = BlockadeShift(atom,nnc6[kk],0,0.5,0.5)
  junk1, junk2, junk3, c6D3_2[kk] = BlockadeShift(atom,nnc6[kk],2,1.5,1.5)
  junk1, junk2, junk3, c6D5_2[kk] = BlockadeShift(atom,nnc6[kk],2,2.5,2.5)
end

save("Data/C6data.jld","nnc6",nnc6,"c6S1_2",c6S1_2,"c6D3_2",c6D3_2,"c6D5_2",c6D5_2)
