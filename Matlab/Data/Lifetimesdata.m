%Getting lifetimes of 87Rb from Branden et al. (2009)

close all
clear
clc

sexpt = [15.2 14.9 18.9 17.9 18.1 20.1 23.8 24.4 25.7 26.8 28.6 34.3...
    34.2 38.1 40.9 40.1 47 50.8];
sexpterr = [0.7 0.8 0.8 0.8 0.5 0.8 1.7 1.0 1.4 1.0 1.1 1.5 0.9 2.1 1.2...
    1.6 1.7 2.4];
pexpt = [33.2 33.7 39.8 43.5 45.1 45.9 52.6 58.7 60.3 60.2 64.2];
pexpterr = [1.4 1.5 1.3 1.9 2.0 1.3 1.3 3.0 1.7 2.8 2.6];
dexpt = [13.5 15.8 19.5 19.9 20.5 24.1 25.4 28.6 31.5 32.3 36.4 35.0...
    36.9 38.7 38.9 42];
dexpterr = [0.7 0.5 1.5 1.2 2.3 0.6 1.8 1.6 0.8 1.2 1.4 1.0 1.8 2.4 3.0...
    3.2];

save('Lifetimes.mat','sexpt','sexpterr','pexpt','pexpterr','dexpt',...
'dexpterr');