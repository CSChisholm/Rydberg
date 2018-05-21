%Numbers for 87Rb
%24/07/2017

%Low level parameters
ionlim = 4.1771270; %Rb groundstate energy in eV
D1trans = 1.56051536; %D1 line transition energy in eV
D2trans = 1.58999085; %D2 line transition energy in eV
n4l2j5_2en = 2.40116159;
n4l2j3_2en = 2.40121692;
n6l0j1_2en = 2.49759249;
n6l1j1_2en = 2.94203794;
n6l1j3_2en = 2.95165365;
n5l0j1_2en = 0; %5S_{1/2} energy in eV
%Other parameters
ZZ = 37; %Nucler charge of Rubidium
spin = 1/2; %Electron spin
%Quantum defects from Li et al. (2003) and Han et al. (2006)
if ll==0
    delta_0 = 3.1311804;
    delta_2 = 0.1784;
elseif ll==1&&jj==1/2
    delta_0 = 2.65448849;
    delta_2 = 0.2900;
elseif ll==1&&jj==3/2
    delta_0 = 2.6416737;
    delta_2 = 0.2950;
elseif ll==2&&jj==3/2
    delta_0 = 1.34809171;
    delta_2 = -0.60286;
elseif ll==2&&jj==5/2
    delta_0 = 1.34646572;
    delta_2 = -0.59600;
elseif ll==3&&jj==5/2
    delta_0 = 0.0165192;
    delta_2 = -0.085;
elseif ll==3&&jj==7/2
    delta_0 = 0.0165437;
    delta_2 = -0.086;
else
    delta_0 = 0;
    delta_2 = 0;
end
%Approximation for quantum efect assuming n > 20
delta_nlj = delta_0 + delta_2/(nn-delta_0)^2;
%Model parameters from Marinescu et al. (1994)
Ry = 1;
alpha_c = 9.0760;
if ll==0
    a_1 = 3.69628474;
    a_2 = 1.64915255;
    a_3 = -9.86069196;
    a_4 = 0.19579987;
    r_c = 1.66242117;
elseif ll==1
    a_1 = 4.44088978;
    a_2 = 1.92828831;
    a_3 = -16.79597770;
    a_4 = -0.81633314;
    r_c = 1.50195124;
elseif 11==2
    a_1 = 3.78717363;
    a_2 = 1.57027864;
    a_3 = -11.65588970;
    a_4 = 0.52942835;
    r_c = 4.86851938;
else
    a_1 = 2.39848933;
    a_2 = 1.76810544;
    a_3 = -12.07106780;
    a_4 = 0.77256589;
    r_c = 4.79831327;
end

%Get energy eigen value                         
if nn==5&&ll==0
    eneigval = -ionlim/atomenergy;
elseif nn==5&&ll==1&&jj==1/2
    eneigval = (D1trans-ionlim)/atomenergy;
elseif nn==5&&ll==1&&jj==3/2
    eneigval = (D2trans-ionlim)/atomenergy;
elseif nn==4&&ll==2&&jj==5/2
    eneigval = (n4l2j5_2en-ionlim)/atomenergy;
elseif nn==4&&ll==2&&jj==3/2
    eneigval = (n4l2j3_2en-ionlim)/atomenergy;
elseif nn==6&&ll==0
    eneigval = (n6l0j1_2en-ionlim)/atomenergy;
elseif nn==6&&ll==1&&jj==1/2
    eneigval = (n6l1j1_2en-ionlim)/atomenergy;
elseif nn==6&&ll==1&&jj==3/2
    eneigval = (n6l1j3_2en-ionlim)/atomenergy;
else
   eneigval = -Ry/(2*(nn-delta_nlj)^2);
end