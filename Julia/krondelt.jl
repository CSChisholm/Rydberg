#Kronecker delta function
#23/07/2017

function krondelt(ii,jj)
  if ii==jj
    kronecker = 1
  else
    kronecker = 0
  end
  return kronecker
end
