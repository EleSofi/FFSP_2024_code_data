function [V,w] = jump(V_,w_,y_,y,C,S_u,S_o,S_bis,mu,mtot,W_star,MAK,network)
%UNTITLED14 Summary of this function goes here
%   Detailed explanation goes here
k=0;
O_k=zeros(1,1);
for j=mu+1:mtot
        if S_o(:,j) == y-y_
           k=k+1;
           if k == 1
               O_k(k)=j;
           else
               O_k=[O_k,j];
           end
        end
end

if MAK == 1 && network == 1
    [W] = propensity_bimolecular(V_,y_,C,S_bis);
elseif MAK == 1 && network == 0
    W=propensity(V_,y_,C,S_bis);
else
    W=zeros(1,length(W_star));
       for k=1:length(W_star)
           W(k)=W_star{k}(V_,y_);
       end
       
end
lambda=zeros(1,length(O_k));
for i=1:length(O_k)
    lambda(i)=W(O_k(i));
end
if length(O_k) > 1
u=rand;
Index=1;
while (Index < u*length(O_k)) && (Index < length(O_k))
    Index=Index+1;
end
else
    Index=1;
end
V=V_+S_u(:,O_k(Index));
w=w_*W(O_k(Index));

end

