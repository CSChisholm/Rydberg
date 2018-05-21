#Function for evaluating radial matrix elements
#23/07/2017

function radiel(atom,nn1,ll1,jj1,nn2,ll2,jj2)

  #Compute radial wavefunctions
    radial1, rscale1 = numerovfunc(atom,nn1,ll1,jj1)

    radial2, rscale2 = numerovfunc(atom,nn2,ll2,jj2)

  #Some preparation
  if nn1>=nn2
    Yk1 = radial1
    Yk2 = radial2
    rscale = rscale1
  else
    Yk1 = radial2
    Yk2 = radial1
    rscale = rscale2
  end

  #Resize smaller solution vector by attaching zeros to the end such that the two solution vectors are the same length
  if length(Yk1)≠length(Yk2)
    Yk2conc = append!(Yk2,zeros(length(Yk1)-length(Yk2)))
  else
    Yk2conc = Yk2
  end

  #Solve the matrix elements using method adapted from Zimmerman et al. (1979)
  deltar = rscale[2:end]-rscale[1:(end-1)]
  numervec = Yk1[2:end].*Yk2conc[2:end].*(rscale[2:end].^2).*deltar
  return matrixelement = sum(numervec)
end
