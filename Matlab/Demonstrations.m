%Script to make a series of demonstration plots
%24/07/2017

close all
clear
clc

%Input information
atom = '87Rb';
nn = 50;
ll = 0;
jj = 0.5;
mj = 0.5;

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

%Plot wave function
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
ylabel('$|rR(r)|.^2\, (a_0^{-1})$','interpreter','latex');
title([atom,' radial wavefunction |n,l,j\rangle = |',num2str(nn),',',...
    num2str(ll),',',num2str(jj),'\rangle']);
set(gcf,'Color','w');

%Compare ground state to self-consistent field method
load('Data/SCFcalc.mat');
[normY_solg,rrg] = numerovfunc(atom,5,0,0.5);

figure
plot(plotscale2, sqrtrR, 'o');
hold on
plot(sqrt(rrg), -normY_solg)
legend('Callaway', 'Numerov');
xlabel('$(r/a_0)^{1/2}$','interpreter','latex');
ylabel('$r^{1/2}R(r)\, (a_0^{-1})$','interpreter','latex');
title('Rubidium Groundstate Radial Wavefunction');
set(gcf,'Color','w');

%Lifetimes
load('Data/Lifetimes.mat');

%S_1/2
nns = 28:45;
calcs = zeros(size(nns));
for kk = 1:length(nns)
    calcs(kk) = Radiative_Lifetimes(atom,nns(kk),0,0.5)*1e6;
end

%P_3/2
nnp = 34:44;
calcp = zeros(size(nnp));
for kk = 1:length(nnp)
    calcp(kk) = Radiative_Lifetimes(atom,nnp(kk),1,1.5)*1e6;
end

%D_5/2
nnd = 29:44;
calcd = zeros(size(nnd));
for kk = 1:length(nnd)
    calcd(kk) = Radiative_Lifetimes(atom,nnd(kk),2,2.5)*1e6;
end

figure
plot(nns, calcs, 'x');
hold on
errorbar(nns, sexpt, sexpterr, 'o');
title('Plot of calculated and experimental lifetimes for nS_{1/2} states');
ylabel('$\tau\, (\mu\mathrm{s})$','interpreter','latex');
xlabel('$n$','interpreter','latex');
legend('Calculated', 'Experimental', 'Location', 'SouthEast');
hold off
set(gcf,'Color','w');

figure
plot(nnp, calcp, 'x');
hold on
errorbar(nnp, pexpt, pexpterr, 'o');
title('Plot of calculated and experimental lifetimes for nP_{3/2} states');
ylabel('$\tau\, (\mu\mathrm{s})$','interpreter','latex');
xlabel('$n$','interpreter','latex');
legend('Calculated', 'Experimental', 'Location', 'SouthEast');
hold off
set(gcf,'Color','w');

figure
plot(nnd, calcd, 'x');
hold on
errorbar(nnd, dexpt, dexpterr, 'o');
title('Plot of calculated and experimental lifetimes for nD_{5/2} states');
ylabel('$\tau\, (\mu\mathrm{s})$','interpreter','latex');
xlabel('$n$','interpreter','latex');
legend('Calculated', 'Experimental', 'Location', 'SouthEast');
hold off
set(gcf,'Color','w');

%Blockade shift

[RRSI,theta,blockadeshiftGHzmesh,C_6val] =...
    BlockadeShift(atom,nn,ll,jj,mj);

stringlookup = 'SPD';
[numer, denom] = rat(jj);
[numer1, denom1] = rat(mj);

figure
imagesc(RRSI*1e6,theta,blockadeshiftGHzmesh*1e3);
xlabel('$R\, (\mu\mathrm{m})$','interpreter','latex');
ylabel('$\theta\, (\mathrm{radians})$','interpreter','latex');
zlabel('$\Delta\!W\, (\mathrm{MHz})$','interpreter','latex');
set(gcf,'Color','w');
axis tight;
colorbar
title(['Calculated Rydberg blockade shift for |', num2str(nn),...
    stringlookup(ll+1), '_{', num2str(numer), '/', num2str(denom),...
    '}, m_j = ', num2str(numer1), '/', num2str(denom1),...
    '\rangle state of Rb']);

%C_6
load('Data/C6data.mat','nnc6','c6S1_2','c6D3_2','c6D5_2');
%Convert from GHz/m^6 to GHz/um^6
c6S1_2 = c6S1_2*1e36;
c6D3_2 = c6D3_2*1e36;
c6D5_2 = c6D5_2*1e36;

figure;
semilogy(nnc6,abs(c6S1_2),'.-');
xlabel('$n$','interpreter','latex');
ylabel('$|C_6|\, (\mathrm{GHz}\cdot\mu\mathrm{m}^6)$','interpreter',...
    'latex');
title('Calculated C_6 for nS_{1/2}');
xlim([nnc6(1),nnc6(end)]);
set(gca,'Yscale','log');
set(gcf,'Color','w');

figure;
semilogy(nnc6,abs(c6D3_2),'.-');
xlabel('$n$','interpreter','latex');
ylabel('$|C_6|\, (\mathrm{GHz}\cdot\mu\mathrm{m}^6)$','interpreter',...
    'latex');
title('Calculated C_6 for nD_{3/2}');
xlim([nnc6(1),nnc6(end)]);
set(gca,'Yscale','log');
set(gcf,'Color','w');

figure;
semilogy(nnc6,abs(c6D5_2),'.-');
xlabel('$n$','interpreter','latex');
ylabel('$|C_6|\, (\mathrm{GHz}\cdot\mu\mathrm{m}^6)$','interpreter',...
    'latex');
title('Calculated C_6 for nD_{5/2}');
xlim([nnc6(1),nnc6(end)]);
set(gca,'Yscale','log');
set(gcf,'Color','w');
