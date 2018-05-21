#Wigner 6j-symbol calculator.
#Translated and modified from function by Amita B Deb, Clarendon Lab. 2007.

# Calculates { j1, j2 ,j3}  using Racah formula. See: Sobelman: Atomic Spectra and Radiative Transitions.
#              J1  J2  J3

function Wigner6j(j1,j2,j3,J1,J2,J3)

  check = Wigner6jcheck(j1,j2,j3,J1,J2,J3)

  if check==1
    #Finding triangular coefficients
    tri1 = triangle_coeff(j1,j2,j3)
    tri2 = triangle_coeff(j1,J2,J3)
    tri3 = triangle_coeff(J1,j2,J3)
    tri4 = triangle_coeff(J1,J2,j3)

    #Finding the range of summation in Racah formula
    a = zeros(4)
    a[1] = j1 + j2 + j3
    a[2] = j1 + J2 + J3
    a[3] = J1 + j2 + J3
    a[4] = J1 + J2 + j3

    rangei = maximum(a)

    k = zeros(3)
    k[1] = j1 + j2 + J1 + J2
    k[2] = j2 + j3 + J2 + J3
    k[3] = j3 + j1 + J3 + J1

    rangef = minimum(k)

    Wignerhold = 0

    for tt = rangei:rangef
      Wignerhold = Wignerhold+ (tri1*tri2*tri3*tri4)^(0.5)*((-1)^tt)*factorial(tt+1)/fung(tt,j1,j2,j3,J1,J2,J3)
    end
    return Wigner = Wignerhold
  else
    return Wigner = 0
  end
end
