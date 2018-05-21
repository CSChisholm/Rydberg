%Data from Callaway and Morgan

close all
clear
clc

rr = [0:0.005:0.02, 0.03:0.01:0.1, 0.12:0.02:0.30, 0.35:0.05:0.6,...
    0.7:0.1:1.2, 1.4:0.2:4.4, 4.8, 5.2, 5.6:0.4:8.4, 9.2:0.8:13.2,...
    14.8:1.6:26];

R_a = [0, 0.02055, 0.03362, 0.04069, 0.043, 0.03745, 0.02353, 0.00576,...
    -0.01281, -0.03018, -0.04517, -0.05711, -0.06573, -0.07309,...
    -0.06882, -0.05565, -0.03648, -0.01396, 0.00967, 0.03269, 0.05386,...
    0.07231, 0.08754, 0.11040, 0.11257, 0.09803, 0.07170, 0.03828,...
    0.00167, -0.06995, -0.12867, -0.1691, -0.1907, -0.1950, -0.1848,...
    -0.1333, -0.0580, 0.0252, 0.1064, 0.1807, 0.2464, 0.3032, 0.3510,...
    0.3902, 0.4215, 0.4456, 0.4632, 0.4751, 0.4819, 0.4844, 0.4831,...
    0.4713, 0.4503, 0.4233, 0.3927, 0.3604, 0.3276, 0.2953, 0.2643,...
    0.2350, 0.2078, 0.1601, 0.1213, 0.0907, 0.0669, 0.0489, 0.0354,...
    0.0180, 0.0089, 0.0044, 0.0021, 0.001, 0.0005, 0.0002, 0.0001];

%Tidy up for plotting

truncrr = rr(2:length(R_a));
plotscale2 = sqrt(truncrr);
truncR_a = R_a(2:length(R_a));
sqrtrR = truncR_a./plotscale2;

save('SCFcalc.mat','plotscale2','sqrtrR');