#Function for getting solution to radial Schrodinger equation using numerov algorithm
#23/07/2017

function numerovfunc(atom,nn,ll,jj)

#Get parameters
ZZ, spin, eneigval, params = GetAtomParams(atom,nn,ll,jj)
alpha_c = params[1]
a_1 = params[2]
a_2 = params[3]
a_3 = params[4]
a_4 = params[5]
r_c = params[6]

#Create vector to integrate over
hh = 0.01 #Choose this to be small so that O(hh^6) may be ignored
h2 = hh^2 #For efficiency
r_i = hh/2 #Define starting point
r_o = 2*nn*(nn+15)  #Define end point
x_i = log(r_i) #Change of variables from Callagher, 2005
x_o = log(r_o)
xx = collect(x_i:hh:x_o)
rr = exp(xx) #For efficiency

#Set up Schrodinger equation from Gallagher (2005), note Pritchard (2012)

LdotS = (jj*(jj+1)-ll*(ll+1)-spin*(spin+1))/2 #Spin-orbit coupling

spinorbitpotential = finestrucconst^2*LdotS./(2*rr.^3) #Fine structure splitting

radialcharge = 1 + (ZZ-1)*exp(-a_1*rr) - rr.*(a_3+a_4*rr).*exp(-a_2*rr) #Effective nuclear charge

coulombpotential = -radialcharge./rr - (alpha_c./(2*(rr.^4))).*(1-exp(-((rr/r_c).^6))) #Corrected Coulomb potential

totalpotential = spinorbitpotential+coulombpotential #total potential

cenfugterm = (ll+1/2)^2 #centrifugal term from Schrodinger equation

#Apply Numerov algorithm
G_x = 2*exp(2*xx).*(totalpotential-eneigval) + cenfugterm #Coeficient in differential equation

T_x = h2*G_x/12 #For efficiency

Ysoln_vec = zeros(size(xx)) #Generate a placeholder vector for solutions

Ysoln_vec[end] = -1e-10 #The function must go to zero at infinity

Ysoln_vec[end-1] = -2e-10

#To perform the iteration with convenient indexing the vectors xx, G_x and soln_vec are flipped

fYsoln = flipdim(Ysoln_vec,1)
fT_x = flipdim(T_x,1)

for ii = 3:length(xx)
        fYsoln[ii] = ((2+10*fT_x[ii-1])*fYsoln[ii-1] -  (1-fT_x[ii-2])*fYsoln[ii-2])/(1-fT_x[ii])
end

Ysoln_vec = flipdim(fYsoln,1) #Return solutions to proper ordering

#Normalise (method adapted from Gallagher, 2005)
normconstvec = zeros(length(rr)-1)
for mm = 2:length(rr)
    deltar = rr[mm]-rr[mm-1]
    normconstvec[mm-1] = (Ysoln_vec[mm]^2)*rr[mm]*deltar
end

normconst = sqrt(sum(normconstvec))
normY_sol = Ysoln_vec/normconst

#Saving and loading seems to hurt performance
#filename = filenamemake(atom,nn,ll,jj)

#save(filename,"rr",rr,"normY_sol",normY_sol)

return normY_sol, rr
end
