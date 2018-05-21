%Function for getting atom specfifc parameters for Rydberg code. Based on
%Matlab code written in summer of 2015 updated 24/07/2017.
%The intention is not that this be a complete database but rather that more
%parameters be added as they are needed

function [ZZ,spin,eneigval,params] = GetAtomParams(atom,nn,ll,jj) %#ok<STOUT,INUSD>
SIunits;
if strcmp(atom,'87Rb')==1
    Rb87numbers;
end
params = [alpha_c,a_1,a_2,a_3,a_4,r_c];