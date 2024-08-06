function [V,w] = Continuous_Evolution(V0,w0,t0,tf,y,C,S_u,mu,mtot,S_bis,W_star,MAK,network)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
V=V0;
w=w0;
t=t0;
i=1;
lambda=zeros(1,mtot);
while t(i) < tf
    if MAK == 1 && network == 1
    [W] = propensity_bimolecular(V(:,i),y,C,S_bis);
    elseif MAK ==1 && network == 0
        [W] = propensity_bimolecular(V(:,i),y,C,S_bis);
    else
        W=zeros(1,length(W_star));
        for k=1:length(W_star)
           W(k)=W_star{k}(V(:,i),y);
        end
    end
    
    
    
    for j=1:mu
        lambda(j)=W(j);
    end
    for j=mu+1:mtot
        lambda(j)=W(j);
    end
    lambda_u=sum(lambda(1:mu));
    lambda_o=sum(lambda(mu+1:end));
    tau=Exponential(lambda_u);
    if t(i)+tau < tf 
        i_next=next_reaction(W(1:mu));
        V=[V,V(:,i)+S_u(:,i_next)];
    end
    tn=min(t(i)+tau,tf);
    w(i+1)=w(i)*exp(-(tn-t(i))*lambda_o);
    i=i+1;
    t(i)=tn;
        
end
V=V(:,2:end);
w=w(2:end);

end

