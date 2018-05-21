#Calculating the denominator in Racah Formula

function fung(tt,j1,j2,j3,J1,J2,J3)

  return r = factorial(tt-j1-j2-j3)*factorial(tt-j1-J2-J3)*factorial(tt-J1-j2-J3)*factorial(tt-J1-J2-j3)*factorial(j1+j2+J1+J2-tt)*
   factorial(j2+j3+J2+J3-tt)*factorial(j3+j1+J3+J1-tt)
 end
