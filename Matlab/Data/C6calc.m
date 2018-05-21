nnc6 = 35:82;

c6S1_2 = zeros(size(nnc6));
c6D3_2 = zeros(size(nnc6));
c6D5_2 = zeros(size(nnc6));

for kk = 1:length(nnc6)
[~,~,~,c6S1_2(kk)] = BlockadeShift(atom,nnc6(kk),0,0.5,0.5);
[~,~,~,c6D3_2(kk)] = BlockadeShift(atom,nnc6(kk),2,1.5,1.5);
[~,~,~,c6D5_2(kk)] = BlockadeShift(atom,nnc6(kk),2,2.5,2.5);
end

save('C6data.mat','nnc6','c6S1_2','c6D3_2','c6D5_2');