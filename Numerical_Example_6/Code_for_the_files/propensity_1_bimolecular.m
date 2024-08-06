function [W] = propensity_1_bimolecular(X,Y,C,S_bis,index)
n1=length(X);
n2=length(Y);

W=C(index);
for j=1:n1
    if S_bis(j,index) == 1 
        W=W*X(j);
    end
    if S_bis(j,index) == 2
        W=W*X(j)*((X(j)-1)/2);
    end
end

for j=1:n2 
    if S_bis(n1+j,index) == 1
        W=W*Y(j);
    end
    if S_bis(n1+j,index) == 2
        W=W*Y(j)*((Y(j)-1)/2);
    end
end
    
end

