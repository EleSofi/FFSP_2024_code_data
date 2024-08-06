function [W] = propensity_bimolecular(X,Y,C,S_bis)
%function which creates mass action kinetics functions as many as
%the columns in S
%S is a matrix with rows equal to the number of species, 
%and columns equal to the reactions
n1=length(X);
n2=length(Y);
Shape_S=size(S_bis);
W=zeros(Shape_S(2),1);
for i=1:length(W)
    W(i)=C(i);
end



for i=1:length(W)
    
    for j=1:n1
        if S_bis(j,i) == 1
            W(i)=W(i)*X(j);
        end
        if S_bis(j,i) == 2
            W(i)=W(i)*((X(j)-1)/2)*X(j);
        end     
    end
    for k=1:n2
        if S_bis(n1+k,i) == 1
            W(i)=W(i)*Y(k);
        end
        if S_bis(n1+k,i) == 2
            W(i)=W(i)*Y(k)*((Y(k)-1)/2);
        end    
    end
    
end
        
end
