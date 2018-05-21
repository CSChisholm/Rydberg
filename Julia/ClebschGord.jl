#Function to calculate Clebsch-Gordan coeficients
#22/01/2016

function ClebschGord(jj1,jj2,JJ,mm1,mm2,MM)
return CG =
    (-1)^(jj1-jj2+MM)*sqrt(2*JJ+1)*Wigner3j(jj1,jj2,JJ,mm1,mm2,-MM)
end
