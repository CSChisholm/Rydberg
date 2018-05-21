#Function for calculating Rydberg blockade shifts
#24/07/2017

function BlockadeShift(atom,nn,ll,jj,mj)

  pceV = 2*π*rpceV #Planck's constant in eV
  ZZ, spin, eneigval, params = GetAtomParams(atom,nn,ll,jj)

  #Create ranges for θ and R (R in atomic units)
  θ = collect(0:0.01:(π/2))
  RRSI = (collect(4:0.01:10))*1e-6
  RR = RRSI/bohrrad

  #Set defect threshold
  defectmax = 100e9 #Maximum energy defect in Hz
  entol = defectmax*pceV/atomenergy #Converted to atomic units
  ntol = 4 #Maximum change in n

  #Determine single particle state energy
  singen = eneigval

  #Create set of allowable interacting state quantum numbers
  nvec = rectpulse((nn-ntol):(nn+ntol),ll+2)
  lcands = collect((ll-1):2:(ll+1))
  lcands = filter(lcands->lcands≥0,lcands)
  if lcands[1]==0
    smalllvec = [0, rectpulse(lcands[2:end],2)]
    jmaker = [0.5, repmat([-0.5, 0.5],1,length(lcands[2:end]))]
  else
    smalllvec = rectpulse(lcands,2)
    jmaker = repmat([-0.5, 0.5],1,length(lcands))
  end
  jmaker = reshape(jmaker,length(smalllvec),1)
  smalljvec = smalllvec+jmaker
  lvec = repmat(smalllvec,1,Int(length(nvec)/length(smalllvec)))
  jvec = repmat(smalljvec,1,Int(length(nvec)/length(smalljvec)))
  lvec = reshape(lvec,length(nvec),1)
  jvec = reshape(jvec,length(nvec),1)
  truncspace1 = hcat(nvec, lvec, jvec)'
  pairs1 = combvechack(truncspace1,truncspace1)
  Sp1 = size(pairs1)

  #Check which of these states satisfy the infinite separation energy defect condition and selection rules for j

  pindex = zeros(Sp1[2])
  infdefects = zeros(Sp1[2])

  for kk = 1:Sp1[2]
    energy1 = envalfunc(atom,pairs1[1,kk],pairs1[2,kk],pairs1[3,kk])
    energy2 = envalfunc(atom,pairs1[4,kk],pairs1[5,kk],pairs1[6,kk])
    infdefects[kk] = energy1 + energy2 - 2*singen #Energy defect
    jvec = collect((jj-1):(jj+1))
    j1check = length(find(jvec -> jvec==pairs1[3,kk],jvec))
    j2check = length(find(jvec -> jvec==pairs1[6,kk],jvec))
    if abs(infdefects[kk])≤entol&&j1check==1&&j2check==1
      pindex[kk] = 1
    end
  end
  qvals = find(pindex -> pindex==1,pindex)
  pairs2 = vcat(pairs1[:,qvals],infdefects[qvals]')
  Lp2 = Int(length(infdefects[qvals]))

  matel2part = zeros(length(θ),Lp2) #Array to store matrix elements in

  #Call function to calculate matrix elements
  for kk = 1:Lp2
    matel2part[:,kk] = matrixel2p(atom,nn,ll,jj,mj,pairs2[1,kk],pairs2[2,kk],pairs2[3,kk],pairs2[4,kk],pairs2[5,kk],pairs2[6,kk],pairs2[7,kk],θ)
  end

  #Compute blockade shift in atomic units
  summation = -sum(matel2part,2)
  maxs = maximum(abs(summation))
  kindex = find(summation -> abs(summation)==maxs,summation) #Finds index of maximum blockade shift

  blockadeshiftau = summation[kindex]./(RR.^6)

  RRmesh = repmat(RR,1,length(summation))'
  summationmesh = repmat(summation,1,length(RR))
  blockadeshiftmesh = summationmesh./(RRmesh.^6)

  #Convert to GHz
  encon = (atomenergy/pceV)*1e-9 #Factor from atomic energy units to GHz
  blockadeshiftGHz = blockadeshiftau*encon
  blockadeshiftmeshGHz = blockadeshiftmesh*encon

  #Curve fitting
  C_6val= blockadeshiftGHz[end]*(RRSI[end]^6) #GHz/m⁶

  return RRSI, θ, blockadeshiftmeshGHz, C_6val
end
