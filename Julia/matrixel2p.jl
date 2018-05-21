#Function for calculating two particle matrix elements required for second order perturbation theory evalutation of the blockaed shift
#25/07/2017

function matrixel2p(atom,nni,lli,jji,mji,nn1,ll1,jj1,nn2,ll2,jj2,defect,theta)

  #Predefinitions for efficiency
  sint = sin(theta)
  cost = cos(theta)
  sin2t = sint.^2
  cos2t = cost.^2

  #First calculate angular factors as vectors [minus,zero,plus] for each single particle interaction stae

  angvec1 = angcalc(lli,jji,mji,ll1,jj1)

  angvec2 = angcalc(lli,jji,mji,ll2,jj2)

  #Now the raidla factors

  matelem1i = radiel(atom,nni,lli,jji,nn1,ll1,jj1)

  matelem2i = radiel(atom,nni,lli,jji,nn2,ll2,jj2)

#Calculate terms in two particle matrix element using Wigner-Eckhart theorem (Pritchard, 2012), (Reinhard et al., 2007)
line1 = angvec1[3]*angvec2[3] + angvec1[1]*angvec2[3] + angvec1[2]*angvec2[2]*(1-3*cos2t)
line2 = -(3/2)*sin2t*(angvec1[3]*angvec2[3] + angvec1[3]*angvec2[1] + angvec1[1]*angvec2[3] + angvec1[1]*angvec2[1])
line3 = -(3/sqrt(2))*sint.*cost*(angvec1[3]*angvec2[2] + angvec1[1]*angvec2[2] + angvec1[2]*angvec2[3] + angvec1[2]*angvec2[1])

matrixelement = matelem1i*matelem2i*(line1+line2+line3)

matrixel2P = (matrixelement.^2)./defect

return matrixel2P
end
