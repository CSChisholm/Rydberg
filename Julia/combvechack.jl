#Function for mimicking Matlab's combvec. Not named combvec because I have seen at least one Julia package with a function called combvec.
#25/07/2017

function combvechack(args...)
  numargs = length(args)
  Mvec = zeros(numargs)
  Nvec = zeros(numargs)
  kk = 1
  for arg in args
    if ndims(arg)==1
      Nvec[kk] = 1
      Mvec[kk] = length(arg)
    else
      sA = size(arg)
      Mvec[kk] = sA[1]
      Nvec[kk] = sA[2]
    end
    kk = kk + 1
  end
  Nprod = prod(Nvec)
  outmatrix = zeros(sum(Mvec),Nprod)
  Mstart = 1
  kk = 1
  for arg in args
    Mend = Int(sum(Mvec[1:kk]))
    val1 = rectpulse(arg',Int(Nprod/prod(Nvec[kk:end])))'
    outmatrix[Mstart:Mend,:] = repmat(val1,1,Int(Nprod/prod(Nvec[1:kk])))
    Mstart = Mend + 1
    kk = kk + 1
  end
  return outmatrix
end
