function [t,t_ob,delta,match,X_tot_prev,X,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK)
%X0 vector of the initial abundance of the different N species.
%S stoichiometry matrix 
%tf final time of the process
%C vector whose components are the rates of the M reactions.

t=zeros(1,1);
t_ob=zeros(1,1);
match=ones(1,1);
X_tot=X0;
Y_tot=X0(n1+1:end);

jumps=0;
delta=zeros(length(X0)-n1,1);

n2=length(X0)-n1;

Shape_S=size(S);
Rows_S=Shape_S(1);
Columns_S=Shape_S(2);
S_o=S((n1+1):Rows_S,:);
j=1;
k=1;

for i=1:Columns_S
    if S_o(:,i) ~= zeros(n2,1)
        O(k)=i;
        k=k+1;
    else
        U(j)=i;
        j=j+1;
    end
end
if j > 1
S_new=zeros(Rows_S,Columns_S);
for i=1:length(U)
S_new(:,i)=S(:,U(i));
end
count=1;
for i=length(U)+1:Columns_S
S_new(:,i)=S(:,O(count));
count=count+1;
end
S=S_new;
S_o=S((n1+1):Rows_S,:);

S_new_bis=zeros(Rows_S,Columns_S);
for i=1:length(U)
S_new_bis(:,i)=S_bis(:,U(i));
end
count=1;
for i=length(U)+1:Columns_S
S_new_bis(:,i)=S_bis(:,O(count));
count=count+1;
end
S_bis=S_new_bis;

C_new=zeros(size(C));
for i=1:length(U)
    C_new(i)=C(U(i));
end
count=1;
for i=length(U)+1:Columns_S
    C_new(i)=C(O(count));
    count=count+1;
end
C=C_new;

[network] = bimolecular(S_bis);

if MAK == 0
    
    for i=1:length(U)
        propensity_new{i}=W_star{U(i)};
    end
    count=1;
    for i=length(U)+1:Columns_S
        propensity_new{i}=W_star{O(count)};
        count=count+1;
    end
    W_star=propensity_new;
end
mu=length(U);
mtot=Columns_S;
U=1:mu;
O=mu+1:mtot;
else
    mu=0;
    mtot=Columns_S;
    U=mu;
    O=mu+1:mtot;
end
    

i=1;
j=1;

while t_ob(j) < tf 
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
        W(k)=W_star{k}(X_tot(1:n1,i),Y_tot(:,j));
        end
        end
       
        
        w0=sum(W);
        t_next=Exponential(w0);
        
        i_next=next_reaction(W);
        

        if ismember(i_next,O) == 1
        i=i+1;
        j=j+1;
        match(j)=i;
        X_tot_prev=X_tot;
        X_tot(:,i)=X_tot(:,i-1)+S(:,i_next);
        Y_tot_prev=Y_tot;
        Y_tot(:,j)=Y_tot(:,j-1)+S_o(:,i_next);
        t(i)=t(i-1)+t_next;
        t_ob(j)=t(i);%t_ob(j-1)+t_next;
        jumps=jumps+1;
        delta(:,jumps)=Y_tot(:,jumps+1)-Y_tot_prev(:,jumps);
        else
        i=i+1;
        X_tot_prev=X_tot;
        X_tot(:,i)=X_tot(:,i-1)+S(:,i_next);
        t(i)=t(i-1)+t_next;
        end
            
end
             

delta=delta(:,1:length(t_ob)-1);
match=match(1:length(match)-1);
t_ob=t_ob(1:length(t_ob)-1);
t=t(1:length(t)-1);
%Y=Y_tot; %Y_tot_prev(:,:);
Y=Y_tot_prev;
X=zeros(n1,length(t_ob));
for i=1:length(t_ob)-1
    X(:,i)=X_tot_prev(1:n1,match(i));
end

end


