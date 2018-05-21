#Function for calculating Rydberg state lifetimes
#23/07/2017

function Radiative_lifetimes(atom,nn,ll,jj)

  temp = 300 #Approximate room temperature in K

  if atom=="87Rb"
    ZZ, spin, eneigval, params, delta_nlj = Rb87numbers(nn,ll,jj)
  end
  dict = load(string("Data/",atom,"levels.jld"))

  #Radiative lifetimes
  if atom=="87Rb"
    loadparam, stateenvec, Rdecayratevec, qvec1, qvec4 = Rb87DecayPaths(ll,jj)
  end

  for qq = qvec1:qvec4
    #Get matrix element
    matrixelement, matrixelementSI = matrix_elements(atom,nn,ll,jj,loadparam[qq,1],loadparam[qq,2],loadparam[qq,3])
    #Compute energies
    enlow = stateenvec[qq]-dict["ionlim"] #Energy of lower state in eV
    eneigvaleV = eneigval*atomenergy #Energy of state |n,l,j⟩ in eV
    #Calculate lifetime excluding black body radiation
    omega = abs((eneigvaleV-enlow)/rpceV) #Angular frequency of transimission
    term1 = (omega^3)/(3*pi*vacpmtvty*rpcJ*(lightc^3)) #splitting terms to make reading easier
    term2 = (2*jj+1)/(2*loadparam[qq,3]+1)

    Rdecayratevec[qq] = term1*term2*matrixelementSI^2 #Radiative decay rate
  end

  Rdecayrate = sum(Rdecayratevec)

  #Account for black body radiation following Beterov et al. (2009)

  neff = nn - delta_nlj #Effective principal quantum number

  if atom=="87Rb"
    BBdecayrate = Rb87blackbody(neff,temp,ll,jj)
  end

  decayrate = Rdecayrate+BBdecayrate #total decay rate
  lifetime = 1/decayrate

  return lifetime
end
