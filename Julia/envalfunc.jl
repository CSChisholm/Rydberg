#Function for quickly determining the energy of a state
#25/07/2017

function envalfunc(atom,nn,ll,jj)
  ZZ, spin, eneigval, params = GetAtomParams(atom,nn,ll,jj)
  return eneigval
end
