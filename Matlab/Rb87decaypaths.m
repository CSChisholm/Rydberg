%Radiative decay pathways for 87Rb
%24/07/2017

function [loadparam,stateenvec,Rdecayratevec,qvec1,qvec4] =...
        Rb87decaypaths(nn,ll,jj) %#ok<INUSL>
SIunits    
Rb87numbers;
loadparam = [5, 1, 0.5; 5, 1, 1.5; 6, 1, 1.5; 6, 1, 0.5; 5, 0, 0.5; 4,...
    2, 1.5; 6, 0, 0.5; 4, 2, 2.5]; %Values to load based on allowed decay
                                   %pathways
stateenvec = [D1trans, D2trans, n6l1j3_2en, n6l1j1_2en, n5l0j1_2en,...
    n4l2j3_2en, n6l0j1_2en, n4l2j5_2en]; %Energy levels
Rdecayratevec = zeros(1,8); %Create a vector to store decay rates in

if ll==0||(ll==2&&jj==1.5)
    qvec1 = 1;
    qvec4 = 4;
elseif ll==1&&jj==1.5
    qvec1 = 5;
    qvec4 = 8;
elseif ll==1&&jj==0.5
    qvec1 = 5;
    qvec4 = 7;
elseif ll==2&&jj==2.5
    qvec1 = 2;
    qvec4 = 3;
end