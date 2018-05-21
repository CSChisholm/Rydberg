%Kronecker delta function
%22/01/2016

function krondelt = krondelt(ii,jj)
if ii==jj
    krondelt = 1;
else
    krondelt = 0;
end