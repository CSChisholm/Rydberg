%Script for plotting blockade shift results
%24/07/2017

close all
clear
clc

atom = '87Rb';
nn = 50;
ll = 0;
jj = 0.5;
mj = 0.5;

[RRSI,theta,blockadeshiftGHzmesh,C_6val] =...
    BlockadeShift(atom,nn,ll,jj,mj);

stringlookup = 'SPD';
[numer, denom] = rat(jj);
[numer1, denom1] = rat(mj);

figure
mesh(RRSI*1e6,theta,blockadeshiftGHzmesh*1e3);
xlabel('$R\, (\mu\mathrm{m})$','interpreter','latex');
ylabel('$\theta\, (\mathrm{radians})$','interpreter','latex');
zlabel('$\Delta\!W\, (\mathrm{MHz})$','interpreter','latex');
set(gcf,'Color','w');
axis tight;
title(['Calculated Rydberg blockade shift for |', num2str(nn),...
    stringlookup(ll+1), '_{', num2str(numer), '/', num2str(denom),...
    '}, m_j = ', num2str(numer1), '/', num2str(denom1),...
    '\rangle state of Rb']);