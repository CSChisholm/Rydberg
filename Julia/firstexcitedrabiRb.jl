#Calculates the power required for a given beam waist to achieve a user defined Rabi frequency for the transition between the first excited state and a Rydberg state |n,l,j⟩.
#Currently only working for 87Rb
#25/07/2017

function firstexcitedrabiRb(nn,ll,jj)

  #Define key constants, one advantage of only writing this script for one atom is that everything can be contained here, a disadvantage is that this section will have to be replaced when functionality is improved

spin = 0.5 #Electron spin
ng = 5 #Principle quantum number corresponding to first excited state
lp = 1 #Orbital angular momentum of P state
J3_2 = 1.5 #Total angular momentum (excluding nuclear coupling and considering D₂ line)
rpcJ = 1.054571596e-34 #Reduced Planck's constant in J
lightc = 299792458 #Speed of light in m⋅s⁻¹
vacpmtvty = 8.854187817e-12 #Vacuum permitivity in F/m
bohrrad = 0.5291772083e-10 #Bohr radius in m
eleccharge = 1.602176462e-19 #Electron charge
atom = "87Rb"

#Define parameters relating to the decomposition of the hyperfine basis
mjvec = [0.5, 1.5]
ClebschGordg = (1/sqrt(2))*[1, -1]
matpart = zeros(2)

#Choose Rabi frequency
RabifreqM = 1 #MHz
Rabifreq = 2*π*RabifreqM*1e6

#Choose polarisation of light
qq = 1

#Choose beamwaist
beamwaist = 30 #um

#Compute ⟨n',l',j',mⱼ'|er|n,l,j,mⱼ⟩
#The radial matrix element
matrixelement = radiel(atom,ng,lp,J3_2,nn,ll,jj)

#The angular part
for kk = 1:2
  m_jprime = mjvec[kk] - qq
  angfac1 = (-1)^(J3_2-mjvec[kk]+spin+jj+1)
  angfac2 = sqrt((2*J3_2+1)*(2*jj+1)*(2*lp+1)*(2*ll+1))
  angfac3 = Wigner6j(J3_2,1,jj,ll,spin,lp)
  angfac4 = Wigner3j(J3_2,1,jj,-mjvec[kk],qq,m_jprime)
  angfac5 = Wigner3j(lp,1,ll,0,0,0)
  angular = angfac1*angfac2*angfac3*angfac4*angfac5
  matpart[kk] = matrixelement*angular*ClebschGordg[kk]
end
matrixfac = sum(matpart)

#Calculate required intensity
term01 = matrixfac^(-2)
term02 = vacpmtvty*(rpcJ^2)*lightc/(2*(eleccharge^2)*(bohrrad^2))
intensity = term01*term02*(Rabifreq^2) #Intensity in W/m^2

#Convert to power
return power = (π/2)*intensity*((beamwaist*1e-6)^2)
end
