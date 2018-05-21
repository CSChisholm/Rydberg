#Function for getting atom specific parameters for Rydberg code. Based on Matlab code written in the summer of 2015 migrated to Julia on 23/07/2017
#The intention is not that this be a complete database but rather that more parameters be added as they are needed

function GetAtomParams(atom,nn,ll,jj)

  #Parameters for Rubidium-8
  if atom=="87Rb"
    ZZ, spin, eneigval, params, delta_nlj = Rb87numbers(nn,ll,jj)
  end

  return ZZ, spin, eneigval, params, delta_nlj
end
