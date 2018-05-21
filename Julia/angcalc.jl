#Function for angular calculation in Blockade shift
#25/07/2017

function angcalc(lli,jji,mji,ll1,jj1)
  spin = 0.5 #Electron spin
  angvec1 = zeros(3)
  for qq = -1:1
    mj1 = mji - qq
    j1vec = collect(-jj1:jj1)
    if length(find(j1vec -> j1vec==mj1,j1vec))==1
      ang11 = (-1)^float(jji-mji+spin+jj1+1)
      ang12 = sqrt((2*jji+1)*(2*jj1+1)*(2*lli+1)*(2*ll1+1))
      ang13 = Wigner6j(jji,1,jj1,ll1,spin,lli)
      ang14 = Wigner3j(jji,1,jj1,-mji,qq,mj1)
      ang15 = Wigner3j(lli,1,ll1,0,0,0)
      angvec1[qq+2] = ang11*ang12*ang13*ang14*ang15 #Shift for Julia indexing
    end
  end
  return angvec1
end
