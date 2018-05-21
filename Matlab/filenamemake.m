%Function for constructing filenames
%24/07/2017

function filename = filenamemake(atom,nn,ll,jj)

%Avoid extra '.' in filename
jstring = strrep(num2str(jj),'.','_');

filename = ['Wavefunctions\',...
    atom,num2str(nn),'n',num2str(ll),'l',jstring,'j.mat'];