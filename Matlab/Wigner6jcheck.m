%Function to check Wigner-6j conditions
%09/02/2016

function check = Wigner6jcheck(j1,j2,j3,J1,J2,J3)

Triad = [j1, j2, j3; j1, J2, J3; J1, j2, J3; J1, J2, j3];
tfvec = zeros(8,1);

%Check integer condition
for kk = 1:4
    tfvec(kk) = floor(sum(Triad(kk,:)))==sum(Triad(kk,:));
end

%Check triangle inequalities
for kk = 1:4
    if abs(Triad(kk,1)-Triad(kk,2))<=Triad(kk,3)&&...
            Triad(kk,3)<=(Triad(kk,1)+Triad(kk,2))
        tfvec(kk+4) = 1;
    end
end

tfsum = sum(tfvec);

check = tfsum==8;