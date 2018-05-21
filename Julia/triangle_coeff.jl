#Function to calculate triangle coefficients for angular momenta
#Based on Matlab code written by Amita B. Deb, Clarendon Lab. 2007.

function triangle_coeff(a,b,c)
  return tri = factorial(a+b-c)*factorial(a-b+c)*factorial(-a+b+c)/factorial(a+b+c+1)
end
