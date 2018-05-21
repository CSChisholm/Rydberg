% Calculating the dem=nominator in Racah Formula

function r = fung(t, j1,j2,j3,J1,J2,J3)

r = factorial(t-j1-j2-j3).*factorial(t-j1-J2-J3).*factorial(t-J1-j2-J3).*factorial(t-J1-J2-j3).*factorial(j1+j2+J1+J2-t).*...
   factorial(j2+j3+J2+J3-t).*factorial(j3+j1+J3+J1-t);