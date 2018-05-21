#Function to emulate Matlab's rectpulse
#25/07/2017

function rectpulse(XX,Nsamp)
  UU = repmat(XX,Nsamp,1)
  sX = size(XX)
  YY = zeros(size(UU))
  for kk = 1:sX[1]
    YY[((kk-1)*Nsamp+1):(kk*Nsamp),:] = UU[kk:sX[1]:end,:]
  end
  return YY
end
