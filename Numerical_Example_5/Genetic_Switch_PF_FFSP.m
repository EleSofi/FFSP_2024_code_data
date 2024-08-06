X0=[1;0;0;0];
%Reaction Rates
C=[1,5,8,2,4,1];
%Stoichiometry Matrix
%S=[-1 1 0 0 0 0; 1 -1 0 0 0 0; 0 0 1 0 -1 0 ; 0 0 0 1 0 -1];
S=[1 -1 0 0 0 0; -1 1 0 0 0 0; 0 0 1 -1 0 0; -1 1 0 0 1 -1];
%S_bis=[1 0 0 0 0 0;0 1 1 0 0 0; 0 0 0 1 1 0; 0 0 0 0 0 1];
S_bis=[0 1 0 0 0 0; 1 0 1 0 0 0 ; 0 0 0 1 1 0; 1 0 0 0 0 1];
N_state=100;
n1=3;
A=zeros(N_state+1,n1);
A(:,n1)=0:N_state;
A(:,1)=1;
B=ones(N_state+1,n1);
B(:,1)=0;
B(:,n1)=0:N_state;
State_Space=vertcat(A,B);

n2=length(X0)-n1;

Shape_State_Space=size(State_Space);
State=Shape_State_Space(1);
[c,index]=intersect(State_Space,X0(1:n1)','rows');
p0=zeros(length(State_Space),1);
p0(index)=1;
N_tot=10000;
MAK=1;
W_star=[];
resample=1;
tf=5;

%% Gillespie

[t,t_ob,delta,match,X_tot_prev,X,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);

%% FFSP2
[T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);

%% particle filter 
[V_tot,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E_pf,Var_pf,SD_pf] = particle_filter_1(t_ob,Y,p0,C,tf,N_tot,State_Space,S,S_bis,resample,n1,W_star,MAK);

%% plot
f=figure;
f.Units='points';
f.OuterPosition=[10 10 2000 450];
freq=4;




subplot(1,3,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlabel('t (s)')
ylabel('Y(t)')
xlim([0 t_ob(end)])
title('Observation Process')
set(gca,'FontSize',20)



subplot(1,3,2)
stairs(t_ob,X(2,:),'b','LineWidth',2)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]],2)',SD_FSP([1 [1:freq:end]],2)','r--','LineWidth',3.5)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],2)',SD_pf([1 [1:freq:end]],1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold off
title('Hidden Process (Activated Gene)')
xlim([0 t_ob(end)])
ylim([0 1.5])
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact Trajectory', 'FFSP' ,'BPF','Location','northwest','Orientation','Horizontal')
set(gca,'FontSize',20)


subplot(1,3,3)
stairs(t_ob,X(3,:),'Color','b','LineWidth',2)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]],3)',SD_FSP([1 [1:freq:end]],3)','r--','LineWidth',3)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],3)',SD_pf([1 [1:freq:end]],3),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold off
ylim([0 7])
xlim([0 t_ob(end)])
title('Hidden Process (mRNA)')
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact Trajectory', 'FFSP' ,'BPF','Location','northwest','Orientation','Horizontal')
set(gca,'FontSize',20) 
