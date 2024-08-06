function [W] = propensity_G(X,C,S_bis)
%function which creates mass action kinetics functions as many as
%the columns in S
%S is a matrix with rows equal to the number of species, 
%and columns equal to the reactions
Shape_S=size(S_bis);
W=zeros(Shape_S(2),1);


for i=1:length(W)
    comp=1;
    for j=1:length(X)
        if X(j) >= S_bis(j,i) && X(j) > 0 && S_bis(j,i) > 0
        %comp=comp*(factorial(X(j))/(factorial(S_bis(j,i))*factorial(X(j)-S_bis(j,i))));
        comp=comp*nchoosek(X(j),S_bis(j,i));
        end
        if X(j) == 0 && S_bis(j,i) > 0
            comp=0;
        end
        
    end
    
    
    W(i)=C(i)*comp;
    
end
        
end

