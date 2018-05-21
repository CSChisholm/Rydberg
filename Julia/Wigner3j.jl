#Wigner 3j-symbol calculator
#17/12/2015

function Wigner3j(aa,bb,cc,alpha,beta,gamma)

# Find triangle coefficient using A. B. Deb (2007) script
triang = triangle_coeff(aa, bb, cc)

#Second term
sqrtterm = sqrt(factorial(aa+alpha)*factorial(aa-alpha)*
    factorial(bb+beta)*factorial(bb-beta)*factorial(cc+gamma)
    *factorial(cc-gamma))

#Finding the range of summation in Racah formula
tmin = zeros(3)
tmin[1] = bb - cc - alpha
tmin[2] = aa + beta - cc
tmin[3] = 0

tmax = zeros(3)
tmax[1] = aa + bb - cc
tmax[2] = aa - alpha
tmax[3] = bb + beta

rangei = maximum(tmin)
rangef = minimum(tmax)

#Sum over the function of t
sumscalar = 0
for tt = rangei:rangef

    sumscalar = sumscalar +
        ((-1)^tt)/tfunction3j(tt,aa,bb,cc,alpha,beta)

end

Wigner = ((-1)^float(aa-bb-gamma))*sqrt(triang)*sqrtterm*sumscalar

return Wigner
end
