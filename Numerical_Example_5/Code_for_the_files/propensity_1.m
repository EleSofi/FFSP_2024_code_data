function [W] = propensity_1(X,Y,C,S_bis,index)
%function which creates mass action kinetics functions as many as
%the columns in S
%S is a matrix with rows equal to the number of species, 
%and columns equal to the reactions
n1=length(X);
n2=length(Y);





    comp=1;
    for j=1:n1
        if X(j) >= S_bis(j,index) && X(j) > 0 && S_bis(j,index) > 0
        %comp=comp*(factorial(X(j))/(factorial(S_bis(j,index))*factorial(X(j)-S_bis(j,index))));
        comp=comp*nchoosek(X(j),S_bis(j,index));
        end
        if X(j) == 0 && S_bis(j,index) > 0
            comp=0;
        end
    end
    for k=1:n2
        if Y(k) >= S_bis(n1+k,index) && Y(k) > 0 && S_bis(n1+k,index) > 0
        %comp=comp*(factorial(Y(k))/(factorial(S_bis(n1+k,index))*factorial(Y(k)-S_bis(n1+k,index)))); 
        comp=comp*nchoosek(Y(k),S_bis(n1+k,index));
        end
        if Y(k) == 0 && S_bis(n1+k,index) > 0
            comp=0;
        end
        
    end
    W=C(index)*comp;

        
end