#function to be called when calculating wigner 3j-symbols
#17/12/2015

function tfunction3j(tt,aa,bb,cc,alpha,beta)
return xoft3j = factorial(tt)*factorial(cc-bb+tt+alpha)*factorial(cc-aa+tt-beta)*factorial(aa+bb-cc-tt)*factorial(aa-tt-alpha)*factorial(bb-tt+beta)
end
