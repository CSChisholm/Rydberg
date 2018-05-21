#Script for saving atomic energy levels
#24/07/2017

using JLD

atom = "87Rb";

D1trans = 1.56051536; #D1 line transition energy in eV
D2trans = 1.58999085; #D2 line transition energy in eV
n4l2j5_2en = 2.40116159;
n4l2j3_2en = 2.40121692;
n6l0j1_2en = 2.49759249;
n6l1j1_2en = 2.94203794;
n6l1j3_2en = 2.95165365;
n5l0j1_2en = 0;
ionlim = 4.1771270; #Rb groundstate energy in eV

save(string(atom,"levels.jld"),"D1trans",D1trans,"D2trans",D2trans,"n4l2j5_2en",n4l2j5_2en,"n4l2j3_2en",n4l2j3_2en,"n6l0j1_2en",n6l0j1_2en,"n6l1j1_2en",n6l1j1_2en,"n6l1j3_2en",n6l1j3_2en,"n5l0j1_2en",n5l0j1_2en,"ionlim",ionlim);
