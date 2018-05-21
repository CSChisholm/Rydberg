#Radiative decay pathways for 87Rb
#23/07/2017

function Rb87DecayPaths(ll1,jj1)

    D1trans = 1.56051536 #D1 line transition energy in eV
    D2trans = 1.58999085 #D2 line transition energy in eV
    n4l2j5_2en = 2.40116159
    n4l2j3_2en = 2.40121692
    n6l0j1_2en = 2.49759249
    n6l1j1_2en = 2.94203794
    n6l1j3_2en = 2.95165365
    n5l0j1_2en = 0
    loadparam = [5 1 0.5; 5 1 1.5; 6 1 1.5; 6 1 0.5; 5 0 0.5; 4 2 1.5; 6 0 0.5; 4 2 2.5] #Values to load based on allowed decay pathways
    stateenvec = [D1trans, D2trans, n6l1j3_2en, n6l1j1_2en, n5l0j1_2en, n4l2j3_2en, n6l0j1_2en, n4l2j5_2en] #Energy levels
    Rdecayratevec = zeros(8) #Create a vector to store decay rates in

    if ll1==0||(ll1==2&&jj1==1.5)
        qvec1 = 1
        qvec4 = 4
    elseif ll1==1&&jj1==1.5
        qvec1 = 5
        qvec4 = 8
    elseif ll1==1&&jj1==0.5
        qvec1 = 5
        qvec4 = 7
    elseif ll1==2&&jj1==2.5
        qvec1 = 2
        qvec4 = 3
    end

    return loadparam, stateenvec, Rdecayratevec, qvec1, qvec4
end
