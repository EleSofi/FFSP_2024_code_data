function [V_tot,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E,Var,SD] = particle_filter_1(t_ob,Y,p0,C,tf,N_tot,State_Space,S,S_bis,resample,n1,W_star,MAK)
%Input Arguments: 
%t_ob observation process jumping times
%Y observation process at jump times
%p0 initial probability distribution of the system
%C Reaction Rates Vector
%N_tot filter sample size
%State_Space State Space of the Hidden Species
%S stoichiometry matrix
%S_bis matrix cointaining the vectors carrying the information 
%about how each species is being consumed in the i-th reaction


%Generate N_tot iid samples from p0 and assign to V1..V_N_tot
[V0] = sampling(p0,State_Space,N_tot);
w=ones(1,N_tot);
for i=1:N_tot
    V_tot{i,1}=V0(:,i);
    w_tot{i,1}=w(i);
    V_jump{i,1}=V0(:,i);
    w_jump{i,1}=w(i);
    %match between the jumps and the continuous evolution
    match_1{i,1}=1;
    match_2{i,1}=1;
end

Shape_S=size(S);
Rows_S=Shape_S(1);
Columns_S=Shape_S(2);
n2=Rows_S-n1;
S_o=S((n1+1):Rows_S,:);
S_u=S(1:n1,:);
[network] = bimolecular(S_bis);
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
else 
    mu=0;
    mtot=Columns_S;
end


    

Shape_State_Space=size(State_Space);
State=Shape_State_Space(1);

t=0;
y=Y(:,1);
resampling=0;

for k=2:length(Y)
    for i=1:N_tot
        %the evolution of the particles is done one by one including the
        %jump and the resampling

           [V_,w_] = Continuous_Evolution_1(V_tot{i,1}(:,end),w_tot{i,1}(end),t,t_ob(k),y,C,S_u,mu,mtot,S_bis,W_star,MAK,network);
           [V,w] = jump(V_(:,end),w_(end),y,Y(:,k),C,S_u,S_o,S_bis,mu,mtot,W_star,MAK,network);
            V_tot{i,1}=[V_tot{i,1},V_,V];
            w_tot{i,1}=[w_tot{i,1},w_,w];
            V_jump{i,k}=V;
            w_jump{i,k}=w;
            match_1{i,1}=[match_1{i,1},length(V_tot{i,1})];
            match_2{i,1}=[match_2{i,1},length(w_tot{i,1})];
     end
        
 
        V_current_jump=zeros(n1,N_tot);
        w_current_jump=zeros(1,N_tot);
        for i=1:N_tot
            V_current_jump(:,i)=V_jump{i,k};
            w_current_jump(i)=w_jump{i,k};
        end
    
    
    %resampling
    if resample == 1 && ((sum(w_current_jump))^2/sum(w_current_jump.^2)) < N_tot %(N_tot*0.95) N_tot
        
        
        [V_tilda,w_tilda] = Offsprings(V_current_jump,w_current_jump);
        resampling=resampling+1;
        for i=1:N_tot
            V_tot{i,1}(:,end)=V_tilda(:,i);
            w_tot{i,1}(end)=w_tilda(i);
            
            V_jump{i,k}=V_tilda(:,i); 
            w_jump{i,k}=w_tilda(i); 
        end
        
    else
        for i=1:N_tot
            w_tot{i,1}(end)=w_current_jump(i)/sum(w_current_jump)*N_tot;
            w_jump{i,k}=w_tot{i,1}(end); 
        end
    end
    t=t_ob(k);
    y=Y(:,k);
end

for i=1:N_tot
    [V_,w_] = Continuous_Evolution_1(V_tot{i,1}(:,end),w_tot{i,1}(end),t_ob(end),tf,Y(:,end),C,S_u,mu,mtot,S_bis,W_star,MAK,network);
    V_tot{i,1}=[V_tot{i,1},V_];
    w_tot{i,1}=[w_tot{i,1},w_];
end
    
pt=zeros(State,length(t_ob));


for i=1:length(t_ob)
    V_current=zeros(n1,N_tot);
    w_current=zeros(1,N_tot);
    for k=1:N_tot
        V_current(:,k)=V_jump{k,i};
        w_current(k)=w_jump{k,i};
    end
    
    for j=1:State
        %[choices,equal]=intersect(V_current',repmat(State_Space(j,:),N_tot,1),'rows');
        lia=ismember(V_current',State_Space(j,:),'rows');
        %w_right=w_current(V_current == State_Space(j,:));
        w_right=w_current(lia);
        if isempty(w_right) ~= 1
            %w_right=w_current(equal);
            pt(j,i)=sum(w_right)/sum(w_current);
        end
    end
end


E=zeros(length(t_ob),n1);
E_square=zeros(length(t_ob),n1);
Var=zeros(length(t_ob),n1);
SD=zeros(length(t_ob),n1);

for i=1:length(t_ob)
    for j=1:n1
    E(i,j)=State_Space(:,j)'*pt(:,i);
    E_square(i,j)=(State_Space(:,j)').^2*pt(:,i);
    Var(i,j)=E_square(i,j)-E(i,j)^2;
    SD(i,j)=sqrt(Var(i,j));
    end
end

end

        
        


