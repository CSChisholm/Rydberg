%Script for plotting Rydberg radial wavefunctions
%24/07/2017

close all
clear
clc

%Input information
atom = '87Rb';
nn = 5;
ll = 0;
jj = 0.5;

%Get some constants
SIunits;
if strcmp(atom,'87Rb')==1
    Rb87numbers;
end

%Calculate wave function
[normY_sol,rr] = numerovfunc(atom,nn,ll,jj);

%Rescale for plotting
plotscale = sqrt(rr);
probamp = (normY_sol.*plotscale).^2;

figure;
plot(plotscale,normY_sol,'LineWidth',2)
if nn > 20
    xlim([sqrt(alpha_c^(1/3)),max(plotscale)]);
end
xlabel('$(r/a_0)^{1/2}$','interpreter','latex');
ylabel('$r^{1/2}R(r)\, (a_0^{-1})$','interpreter','latex');
title([atom,' radial wavefunction |n,l,j\rangle = |',num2str(nn),',',...
    num2str(ll),',',num2str(jj),'\rangle']);
set(gcf,'Color','w');

figure;
plot(rr,probamp,'LineWidth',2)
if nn > 20
    xlim([alpha_c^(1/3),max(rr)]);
end
xlabel('$r/a_0$','interpreter','latex');
ylabel('$|r^{1/2}R(r)|.^2\, (a_0^{-1})$','interpreter','latex');
title([atom,' radial wavefunction |n,l,j\rangle = |',num2str(nn),',',...
    num2str(ll),',',num2str(jj),'\rangle']);
set(gcf,'Color','w');