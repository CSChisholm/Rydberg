%Matrix elements with 5P_{3/2} and laser intensities with a user defined
%Rabi frequency only for 87Rb
%17/12/2015

function power = firstexcitedrabiRb(nn,ll,jj)

%Define key constants
spin = 0.5; %electron spin
ng = 5; %Principle qunatum number
lp = 1; %Azimuthal quantum number for p state
J3_2 = 1.5; %Total angular momentum of nP_{3/2} state
rpcJ = 1.054571596e-34; %Reduced Planck constant in J s
lightc = 2.99792458e8; %Speed of light in m/s
vacpmtvty = 8.854187817e-12; %Vacuum permitivity in F/m
bohrrad = 0.5291772083e-10; %Bohr radius in m
eleccharge = 1.602176462e-19; %electron charge in C
atom = '87Rb';

%Define parameters relating to the decomposition of the hyperfine basis
mjvec = [0.5, 1.5];
ClebschGordg = 1/sqrt(2)*[1, -1];
matpart = [0, 0];

%Choose Rabi frequency
RabifreqM = 1; %MHz
Rabifreq = 2*pi*RabifreqM*1e6;

%Choose polarisation of light
qq = 1;

beamwaist = 30; %um

matrixelement = radiel(atom,ng,lp,J3_2,nn,ll,jj);

 for kk = 1:2   
    %Calculate the <Jm_J|er|J'm_J'> matrix element
    m_jprime = mjvec(kk)-qq;
    angfac1 = (-1)^(J3_2-mjvec(kk)+spin+jj+1);
    angfac2 = sqrt((2*J3_2+1)*(2*jj+1)*(2*lp+1)*(2*ll+1));
    angfac3 = Wigner6j(J3_2,1,jj,ll,spin,lp);
    angfac4 = Wigner3j(J3_2,1,jj,-mjvec(kk),qq,m_jprime);
    angfac5 = Wigner3j(lp,1,ll,0,0,0);
    angular = angfac1*angfac2*angfac3*angfac4*angfac5;
    matpart(kk) = matrixelement*angular*ClebschGordg(kk);
end

matrixfac = sum(matpart);

%Calculate required intensity
term01 = 1/(matrixfac^2);
term02 = vacpmtvty*rpcJ^2*lightc/(2*eleccharge^2*bohrrad^2);
intensity = term01*term02*Rabifreq^2; %Intensity in W/m^2

power = (pi/2)*intensity*(beamwaist*1e-6)^2;