%Function for calculating dipole matrix elements based on results from
%Rydberg_Wavefunctions_Numerov.m
%24/07/2017

function [matrixelement,matrixelementSI] =...
    matrix_elements(atom,nn1,ll1,jj1,nn2,ll2,jj2)

SIunits;

matrixelement = radiel(atom,nn1,ll1,jj1,nn2,ll2,jj2);
matrixelementSI = matrixelement*bohrrad*eleccharge;