function [T,F,jump_times,E,Var,SD,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C,Hidden_Species,delta,S,S_bis,n1,W_star,MAK)
%Function which produces the truncated conditional probability distribution
%solution to a filtering problem for general chemical reaction networks

%code lines needed to sort the reaction channels according to the
%sets U and O: the stoichiometry matrices rows need to be sorted
%with the hidden species first and then the observed ones

Shape_S=size(S);
Rows_S=Shape_S(1);
Columns_S=Shape_S(2);
n2=Rows_S-n1;
S_u=S(1:n1,:);
S_o=S((n1+1):Rows_S,:);
j=1;
k=1;

for i=1:Columns_S   %observable and unobservale reaction channels 
    if S_o(:,i) ~= zeros(n2,1)
        O(k)=i;
        k=k+1;
    else
        U(j)=i;
        j=j+1;
    end
end

[network] = bimolecular(S_bis); %function which classifies a reaction network


if j > 1 %if there are unobserved species
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
S_u=S(1:n1,:);
S_o=S((n1+1):Rows_S,:);
mu=length(U);
mtot=Columns_S;

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

if MAK == 0   %no mass action kinetics
    
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

[New,index_state] = indeces_FFSP(Hidden_Species,mu,S_u); %function which sorts the state space
%according to the unobserved stoiochiometry vectors

%A=@(Y)Evolution_Matrix_FFSP(C,Hidden_Species,Y,S_bis,mu,S_u,MAK,W_star);
U=1:mu;
O=mu+1:mtot;
else
    mu=0;
    mtot=Columns_S;
    U=mu;
    O=mu+1:mtot;
    New=[];
    index_state=[];
    
end

Shape_State=size(Hidden_Species);
Row_State=Shape_State(1);

for i=1:length(t_ob)-1
    
    tspan=[t_ob(i),t_ob(i+1)];
    if i == 2
    F=p;
    T=t_new';
    end
    
    if i == 1
       
    %[t_new,y]=ode45(@(t,y) Evolution_Probability_FFSP(t,y,C,Hidden_Species,Y(:,i),S_bis,mu,S_u,MAK,W_star), tspan, p0);
    [t_new,y]=ode45(@(t,y) Evolution_Probability_FFSP_2(t,y,C,Hidden_Species,Y(:,i),S_bis,mu,MAK,W_star,New,index_state,network), tspan, p0);
    %[t_new,y]=ode45(@(t,y) Evolution_Probability_FFSP_1(t,y,A(Y(:,i))), tspan, p0);
    else
        
    %[t_new,y]=ode45(@(t,y) Evolution_Probability_FFSP(t,y,C,Hidden_Species,Y(:,i),S_bis,mu,S_u,MAK,W_star), tspan, p(:,end));
    %[t_new,y]=ode45(@(t,y) Evolution_Probability_FFSP_1(t,y,A(Y(:,i))), tspan, p(:,end));
    [t_new,y]=ode45(@(t,y) Evolution_Probability_FFSP_2(t,y,C,Hidden_Species,Y(:,i),S_bis,mu,MAK,W_star,New,index_state,network), tspan, p(:,end));
    T=[T t_new(2:end)'];
    end  
    y=y';
    Shape_Sol=size(y);
    Row=Shape_Sol(1);
    Col=Shape_Sol(2);
    
    k=0;
    for j=1:length(O)
        if S_o(:,O(j)) == delta(:,i) %jump of the observation at time tk
            k=k+1;
            O_k(k)=O(j);
        end
    end
    
    
    
    f=zeros(Row_State,length(O_k));
    W_TOT=zeros(Row_State,length(O_k));   
        
    for j=1:length(O_k)
        
        X_new=zeros(size(Hidden_Species));
        for h=1:Row_State
            X_new(h,:)=Hidden_Species(h,:)-S_u(:,O_k(j))';
        end
        
        for m=1:Row_State
            for n=1:Row_State
                if X_new(m,:) == Hidden_Species(n,:)
                f(m,j)=n;
                end
            end  
        end
        
        for a=1:Row_State
            if f(a,j) > 0
                if MAK == 1 && network == 1
                    %[W] = propensity(X_new(a,:),Y(:,i),C,S_bis);
                    W = propensity_1_bimolecular(X_new(a,:),Y(:,i),C,S_bis,O_k(j));
                elseif MAK ==1 && network == 0
                    %W=zeros(1,length(W_star));
                    %for k=1:length(W_star)
                        %W(k)=W_star{k}(X_new(a,:),Y(:,i));
                    %end
                    W=propensity_1(X_new(a,:),Y(:,i),C,S_bis,O_k(j));
                else
                    inputs = num2cell([X_new(a,:),Y(:,i)]);
                    W=W_star{O_k(j)}(inputs{:});
                end
                %W_TOT(a,j)=W(O_k(j));
                W_TOT(a,j)=W;
                
            end
        end
        
    end
    
    p_new=zeros(Row_State,1);
    for w=1:Row_State
    s=0;
    for g=1:length(O_k)
        if f(w,g) > 0
            s=s+W_TOT(w,g)*y(f(w,g),Col);
        end
    end
    p_new(w)=s;
    end
    
    %p_new=p_new/length(O_k);
    Sum=sum(p_new);
    
    if Sum ~= 0
    p_new=p_new/Sum;
    y(:,end)=p_new;
    end
    p=y;
    if i > 1
        F= [ F p(:,2:end) ];
    end
    
end

v=1;
jump_times=zeros(length(t_ob),1);
for i=1:length(t_ob)
    while t_ob(i) ~= T(v)
        v=v+1;
    end
    jump_times(i)=v;
end

E=zeros(length(t_ob),n1);
E_square=zeros(length(t_ob),n1);
Var=zeros(length(t_ob),n1);
SD=zeros(length(t_ob),n1);
Err_jump=zeros(1,length(t_ob));
rho=zeros(Row_State,length(t_ob));

%statistics at jump times
for i=1:length(t_ob)
    for j=1:n1
    E(i,j)=Hidden_Species(:,j)'*F(:,jump_times(i));
    E_square(i,j)=(Hidden_Species(:,j)').^2*F(:,jump_times(i));
    Var(i,j)=E_square(i,j)-E(i,j)^2;
    SD(i,j)=sqrt(Var(i,j));
    end
    rho(:,i)=F(:,jump_times(i));
    Err_jump(i)=1-sum(F(:,jump_times(i)));
end

E_tot=zeros(length(T),n1);
E_square_tot=zeros(length(T),n1);
Var_tot=zeros(length(T),n1);
SD_tot=zeros(length(T),n1);
v=1;
for i=1:length(T)
    if i == jump_times(v)
        v=v+1;
        for j=1:n1
          E_tot(i,j)=Hidden_Species(:,j)'*F(:,i);
          E_square_tot(i,j)=(Hidden_Species(:,j)').^2*F(:,i);
          Var_tot(i,j)=E_square_tot(i,j)-E_tot(i,j)^2;
          SD_tot(i,j)=sqrt(Var_tot(i,j));
        end
    else
    F_temp=F(:,i)/sum(F(:,i));
    for j=1:n1
    E_tot(i,j)=Hidden_Species(:,j)'*F_temp;%F(:,i);
    E_square_tot(i,j)=(Hidden_Species(:,j)').^2*F_temp; %F(:,i);
    Var_tot(i,j)=E_square_tot(i,j)-E_tot(i,j)^2;
    SD_tot(i,j)=sqrt(Var_tot(i,j));
    end
    end
    
end




end