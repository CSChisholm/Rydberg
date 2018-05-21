% Wigner 6j-symbol calculator.
% Written by Amita B Deb, Clarendon Lab. 2007. Modified by Craig S
% Chisholm, Kjaergaard Lab. 2017

% Calculates { j1, j2 ,j3}  using Racah formula. See: Sobelman: Atomic Spectra and Radiative Transitions.
%              J1  J2  J3


function Wigner = Wigner6j(j1,j2,j3,J1,J2,J3)

check = Wigner6jcheck(j1,j2,j3,J1,J2,J3);

if check==1;
% Finding Triangular coefficients
tri1 = triangle_coeff(j1,j2,j3);
tri2 = triangle_coeff(j1,J2,J3);
tri3 = triangle_coeff(J1,j2,J3);
tri4 = triangle_coeff(J1,J2,j3);

% Finding the range of summation in Racah formula.

a(1) = j1 + j2 + j3;
a(2) = j1 + J2 + J3;
a(3) = J1 + j2 + J3;
a(4) = J1 + J2 + j3;

rangei = max(a);

k(1) = j1 + j2 + J1 + J2;
k(2) = j2 + j3 + J2 + J3;
k(3) = j3 + j1 + J3 + J1;

rangef = min(k);

%range = min([j1+j2-j3 j1+J2-J3 J1+j2-J3 J1+J2-j3 j2+j3-j1 J2+J3-j1 j2+J3-J1 J2+j3-J1 j3+j1-j2 J3+j1-J2 J3+J1-j2 j3+J1-J2])





tt = rangei:rangef;
    
    
    
Wigner = (tri1*tri2*tri3*tri4)^(0.5)*((-1).^tt).*factorial(tt+1)./fung(tt,j1,j2,j3,J1,J2,J3);
Wigner = sum(Wigner);
    
else
    Wigner = 0;
end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
%    % Calculateing triange coefficients for Angular momenta.
% 
% function tri = triangle_coeff(a,b,c)
% 
% tri = factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/(factorial(a+b+c+1));
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% % Calculating the dem=nominator in Racah Formula
% 
% function r = fung(t, j1,j2,j3,J1,J2,J3)
% 
% r = factorial(t-j1-j2-j3)*factorial(t-j1-J2-J3)*factorial(t-J1-j2-J3)*factorial(t-J1-J2-j3)*factorial(j1+j2+J1+J2-t)*...
%    factorial(j2+j3+J2+J3-t)*factorial(j3+j1+J3+J1-t);


