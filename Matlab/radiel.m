%function to evaluate radial matrix elements
%22/01/2016

function radiel = radiel(atom,nn1,ll1,jj1,nn2,ll2,jj2)
%first load radial wavefunctions for calculating the radial part of the
%matrix element
if exist(filenamemake(atom,nn1,ll1,jj1),'file')==2
        load(filenamemake(atom,nn1,ll1,jj1), 'normY_sol','rr');
        radial1 = normY_sol;
        rscale1 = rr;
else
    [radial1, rscale1] = numerovfunc(atom,nn1,ll1,jj1);
end

if exist(filenamemake(atom,nn2,ll2,jj2),'file')==2
        load(filenamemake(atom,nn2,ll2,jj2), 'normY_sol','rr');
        radial2 = normY_sol;
        rscale2 = rr;
else
    [radial2, rscale2] = numerovfunc(atom,nn2,ll2,jj2);
end

%Calculate the radial matrix element
if nn1 >= nn2
    Yk1 = radial1;
    Yk2 = radial2;
    rscale = rscale1;
else
    Yk1 = radial2;
    Yk2 = radial1;
    rscale = rscale2;
end

%Resize smaller solution vector by attaching zeros to the end such that the
%two solution vectors are the same length
if length(Yk1)~=length(Yk2)
    szero = zeros(1,length(Yk1)-length(Yk2));
    Yk2conc = [Yk2, szero];
else
    Yk2conc = Yk2;
end

%Solve the matrix elements using method adapted from Zimmerman et al.
%(1979)
deltar = rscale(2:end)-rscale(1:(end-1));
numervec = Yk1(2:end).*Yk2conc(2:end).*(rscale(2:end).^2).*deltar;
radiel = sum(numervec);