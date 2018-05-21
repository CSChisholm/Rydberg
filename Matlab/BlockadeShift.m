%Function for calculating Rydberg blockade shifts.
%Based on code written on 14/01/2016.
%24/07/2017

function [RRSI,theta,blockadeshiftGHzmesh,C_6val] =...
    BlockadeShift(atom,nn,ll,jj,mj)

%Constants
SIunits;
pceV = 2*pi*rpceV; %Planck's constant in eV
if strcmp(atom,'87Rb')==1
    Rb87numbers;
end

%Create ranges for theta and R (R in atomic units)
theta = 0:0.01:pi/2;
RRSI = (4:0.001:10)*1e-6;
RR = RRSI/bohrrad;

%Set defect threshold
defectmax = 100e9; %Maximum energy defect in Hz
entol = defectmax*pceV/atomenergy; %Converted to atomic units
ntol = 4; %Maximum change in n

%Determine single particle state energy
singen = eneigval;

%Create set of allowable interacting state quantum numbers'
nvec = rectpulse((nn-ntol):(nn+ntol),ll+2);
lcands = (ll-1):2:(ll+1);
lcands = lcands(lcands>0|lcands==0);
if lcands(1)==0
    smalllvec = [0,reshape(rectpulse(lcands(2:end),2),1,2*length(lcands(2:end)))];
    jmaker = [0.5,repmat([-0.5,0.5],1,length(lcands(2:end)))];
else
    smalllvec = rectpulse(lcands,2);
    jmaker = repmat([-0.5,0.5],1,length(lcands));
end
smalllvec = reshape(smalllvec,1,numel(smalllvec));
smalljvec = smalllvec-jmaker;
lvec = repmat(smalllvec, 1, length(nvec)/length(smalllvec));
jvec = repmat(smalljvec, 1, length(nvec)/length(smalljvec));
truncspace1 = [nvec; lvec; jvec];
pairs1 = combvec(truncspace1,truncspace1);

%Check which of these states satisfy the infinite separation energy defect
%condition and selection rules for j

pindex = zeros(1,length(pairs1)); %An empty matrix to store allowed pairs in
defects = zeros(1,length(pairs1)); %Energydefects;

for kk = 1:length(pairs1)
    energy1 = envalfunc(atom,pairs1(1,kk),pairs1(2,kk),pairs1(3,kk));
    energy2 = envalfunc(atom,pairs1(4,kk),pairs1(5,kk),pairs1(6,kk));
    defects(kk) = energy1 + energy2 - 2*singen; %Energy defect
    j1check = sum(pairs1(3,kk)==(jj-1):(jj+1));
    j2check = sum(pairs1(6,kk)==(jj-1):(jj+1));
    if abs(defects(kk))<=entol&&j1check==1&&j2check==1
        pindex(kk) = 1;
    end
end
pindex = find(pindex==1);
pairs2 = [pairs1(:,pindex);defects(pindex)];

matel2part = zeros(length(theta),length(pairs2)); %Vector to store matrix
                                                  %elements in

%Call function to calculate matrix elements
for kk = 1:length(pairs2)
    matel2part(:,kk) = matrixel2p(atom,nn,ll,jj,mj,...
        pairs2(1,kk),pairs2(2,kk),pairs2(3,kk),pairs2(4,kk),...
        pairs2(5,kk),pairs2(6,kk),pairs2(7,kk),theta);
end

%Compute the blockade shift in atomic units
summation = -sum(matel2part,2)';
kindex = abs(summation)==max(abs(summation)); %Finds index of maximum
                                              %blockade shift
blockadeshiftau = summation(kindex)./(RR.^6);

[RRmesh,summationmesh] = meshgrid(RR,summation);
blockadeshiftaumesh = summationmesh./(RRmesh.^6);

%Convert units to GHz
encon = (atomenergy/pceV)*1e-9; %Factor from atomic energy units to GHz
blockadeshiftGHz = blockadeshiftau*encon;
blockadeshiftGHzmesh = blockadeshiftaumesh*encon;

%C6
C_6val = blockadeshiftGHz(end)*(RRSI(end))^6; %GHz/m^6