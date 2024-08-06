X0=[5;0];
%Reaction Rates
C=[40;4;1];
%Stoichiometry Matrix
S=[ 0 1 -1; 1 0 0];
S_bis=[ 1 0 1; 0 0 0];
N_state=500;
n1=1;
n2=length(X0)-n1;
[State_Space] = Hidden_State(0:N_state);
Shape_State_Space=size(State_Space);
State=Shape_State_Space(1);
[c,index]=intersect(State_Space,X0(1:n1)','rows');
p0=zeros(length(State_Space),1);
p0(index)=1;
N_tot=10000;
MAK=0;
W_star{1}=@(X,Y)C(1)*(X./(20+X));
W_star{2}=@(X,Y)C(2);
W_star{3}=@(X,Y)C(3)*X;
resample=0;
tf=2;

%% Gillespie trajectories 

[t,t_ob,delta,match,X_tot_prev,X,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);

%% FFSP1
[T,F,jump_times,E,Var,SD,Err_jump,E_tot,SD_tot,rho,eps] = FFSP_2_error(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);

%% FFSP2
[T_1,P,jump_times_1,E_1,Var_1,SD_1,Err_jump_1,E_tot_1,SD_tot_1,Err_tot_1] = FFSP_2_error_new(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);

%% plot 

f=figure;
f.Units='points';
f.OuterPosition=[10 10 1000 500];
freq=1;
subplot(1,2,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlabel('t (s)')
ylabel('Y(t)')
xlim([0 t_ob(end)])
title('Observation Process')
set(gca,'FontSize',20)  


freq1=1;
subplot(1,2,2)
%stairs(t(t<=t_ob(end)),X_1_tot(1,t<=t_ob(end)),'b','LineWidth',2)
stairs(t,X_tot_prev(1,:),'b','LineWidth',2)
hold on
%errorbar(T([1 [1:freq:end]]),E_tot([1 [1:freq:end]],1)',SD_tot([1 [1:freq:end]],1)','r--','LineWidth',5)
errorbar(t_ob([1 [1:freq:end]]),E([1 [1:freq:end]],1)',SD([1 [1:freq:end]],1)','r--','LineWidth',3)
hold on
errorbar(t_ob([1 [1:freq:end]]),E_1([1 [1:freq:end]],1)',SD_1([1 [1:freq:end]],1)','--','Color',[0.4660 0.6740 0.1880],'LineWidth',2)
%errorbar(T_1([1 [1:freq1:end]]),E_tot_1([1 [1:freq1:end]],1)',SD_tot_1([1 [1:freq1:end]],1)','--','Color',[0.4660 0.6740 0.1880],'LineWidth',3)
% hold on
% errorbar(t_ob([1 [1:freq:end]]),E_pf([1 [1:freq:end]],2)',SD_pf([1 [1:freq:end]],1),'--','Color',[0.9290 0.6940 0.1250],'LineWidth',3)
hold off
title('Hidden Process')
xlim([0 t_ob(end)])
ylim([0 10])
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact Trajectory', 'FFSP 1' ,'FFSP 2','Location','northwest')
set(gca,'FontSize',20)  

%%

%FilterDiff=abs(E-E_1);
ProbDiff=zeros(1,length(t_ob));
for i=1:length(t_ob)
    ProbDiff(i)=sum(abs(F(:,jump_times(i))-P(:,jump_times_1(i))));
end
Tab=[t_ob;eps;Err_jump_1;ProbDiff];
csvwrite('/Users/edambrosio/Desktop/PaperFigures/Filter_1_Filter_2_big_par',Tab);

