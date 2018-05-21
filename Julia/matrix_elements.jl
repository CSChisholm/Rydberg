#Function for calculating matrix elements ⟨n₁,l₁,j₁|er|n₂,l₂,j₂⟩ between Rydberg states for a given atom.
#23/07/2017

function matrix_elements(atom,nn1,ll1,jj1,nn2,ll2,jj2)

  matrixelement = radiel(atom,nn1,ll1,jj1,nn2,ll2,jj2)
  matrixelementSI = matrixelement*bohrrad*eleccharge

  return matrixelement, matrixelementSI
end
