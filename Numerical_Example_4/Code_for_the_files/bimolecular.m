function [network] = bimolecular(S_bis)
%This is a function checking if the chemical reaction network is
%bimolecular
network=1;
count=0;
shape_S=size(S_bis);
for i=1:shape_S(2)
    for k=1:shape_S(1)
        if S_bis(k,i) ~= 0
            count=count+1;
        end
    end
    if count > 2 
        network=0;
        return
    end
    count=0;
end
            
end

