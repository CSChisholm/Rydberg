%Function for calculating Rydberg state lifetimes
%24/07/2017

function lifetime = Radiative_Lifetimes(atom,nn,ll,jj)

SIunits;
temp = 300; %Approximation of room temperature in K

if strcmp(atom,'87Rb')==1
    Rb87numbers;
end

%Radiative lifetimes
if strcmp(atom,'87Rb')==1
    [loadparam,stateenvec,Rdecayratevec,qvec1,qvec4] =...
        Rb87decaypaths(nn,ll,jj);
end

for qq = qvec1:qvec4
    %Get matrix element
    [~,matrixelementSI] = matrix_elements(atom,nn,ll,jj,...
        loadparam(qq,1),loadparam(qq,2),loadparam(qq,3));
    %Compute energies
    enlow = stateenvec(qq)-ionlim; %Energy of lower state in eV
    eneigvaleV = eneigval*atomenergy; %Energy of state |nlj> in eV
    
    %Perform calculation
    omega = abs((eneigvaleV-enlow)/rpceV); %Angular frequency of
                                           %transimission
    term1 = omega^3/(3*pi*vacpmtvty*rpcJ*lightc^3); %splitting terms to
                                                    %make reading easier
    term2 = (2*jj+1)/(2*loadparam(qq,3)+1);
    Rdecayratevec(qq) = term1*term2*matrixelementSI^2; %Radiative decay
                                                       %rate
end

Rdecayrate = sum(Rdecayratevec);

%Account for black body radiation Beterov et al. (2009)

neff = nn - delta_nlj; %Effective principal quantum number

if strcmp(atom,'87Rb')==1
    BBdecayrate = Rb87blackbody(neff,temp,ll,jj);
end

decayrate = Rdecayrate+BBdecayrate; %Total decay rate
lifetime = 1/decayrate;