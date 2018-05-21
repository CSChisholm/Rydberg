   % Calculateing triange coefficients for Angular momenta.

function tri = triangle_coeff(a,b,c)

tri = factorial(a+b-c)*factorial(a-b+c)*...
    factorial(-a+b+c)/(factorial(a+b+c+1));