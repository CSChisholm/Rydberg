#Function for constructing filenames
#23/07/2017

function filenamemake(atom,nn,ll,jj)

  #Avoid extra "." in filename
  jstring = replace(string(jj),".","_")

  filename = string("WaveFunctions/",atom,"n",string(nn),"l",string(ll),"j",jstring,".jld")

end
