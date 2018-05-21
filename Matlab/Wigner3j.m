%Wigner 3j-symbol calculator
%17/12/2015

function Wigner = Wigner3j(aa,bb,cc,alpha,beta,gamma)

% Find triangle coefficient using A. B. Deb (2007) script
triang = triangle_coeff(aa, bb, cc);

%Second term
sqrtterm = sqrt(factorial(aa+alpha)*factorial(aa-alpha)*...
    factorial(bb+beta)*factorial(bb-beta)*factorial(cc+gamma)...
    *factorial(cc-gamma));

%Finding the range of summation in Racah formula
tmin(1) = bb - cc - alpha;
tmin(2) = aa + beta - cc;
tmin(3) = 0;

tmax(1) = aa + bb - cc;
tmax(2) = aa - alpha;
tmax(3) = bb + beta;

rangei = max(tmin);
rangef = min(tmax);

%Sum over the function of t

tt = rangei:rangef;

sumscalar = ((-1).^tt)./tfunction3j(tt,aa,bb,cc,alpha,beta,gamma);

Wigner = (-1)^(aa-bb-gamma)*sqrt(triang)*sqrtterm*sum(sumscalar);