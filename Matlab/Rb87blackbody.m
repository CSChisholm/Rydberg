%Function for calculating black body decay rate for 87Rb Rydberg states
%24/07/2017

function BBdecayrate = Rb87blackbody(neff,temp,ll,jj)

%Numbers from Beterov et a;. (2009)
if ll==0
    AA = 0.134;
    BB = 0.251;
    CC = 2.567;
    DD = 4.426;
elseif ll==1&&jj==0.5
    AA = 0.053;
    BB = 0.128;
    CC = 2.183;
    DD = 3.989;
elseif ll==1&&jj==1.5
    AA = 0.046;
    BB = 0.109;
    CC = 2.085;
    DD = 3.901;
elseif ll==2&&jj==1.5
    AA = 0.033;
    BB = 0.084;
    CC = 1.912;
    DD = 3.716;
elseif ll==2&&jj==2.5
    AA = 0.032;
    BB = 0.082;
    CC = 1.898;
    DD = 3.703;
end

term3 = AA/(neff^DD);
term4 = 2.14e10/(exp(315780*BB/((neff^CC)*temp))-1);
BBdecayrate = term3*term4;