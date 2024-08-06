function [t,X] = Gillespie_General(X0,S,S_bis,t,tf,C,W_star,MAK)
%X0 vector of the initial abundance of the different N species.
%S stoichiometry matrix 
%tf final time of the process
%C vector whose components are the rates of the M reactions.



X_tot=X0;






[network] = bimolecular(S_bis);



i=1;

while t(i) < tf 
        %computation of the M propensity functions for
        %the specific M reactions
        %[W] = propensity_G(X_tot(:,i),C,S_bis);
        if MAK == 1 && network == 1 
        [W] = propensity_bimolecular_FSP(X_tot(:,i),C,S_bis);
        elseif MAK == 1 && network == 0
        [W] = propensity_G_FSP(X_tot(:,i),C,S_bis);
        else
        W=zeros(1,length(W_star));
        for k=1:length(W_star)
        W(k)=W_star{k}(X_tot(:,i));
        end
        end
       
        
        w0=sum(W);
        t_next=Exponential(w0);
        
        i_next=next_reaction(W);
        
        
        i=i+1;
        X_tot_prev=X_tot;
        X_tot(:,i)=X_tot(:,i-1)+S(:,i_next);
        t(i)=t(i-1)+t_next;
        
            
end
             


t=t(1:length(t)-1);
t=[t tf];

X=[X_tot_prev X_tot_prev(:,end)];




end