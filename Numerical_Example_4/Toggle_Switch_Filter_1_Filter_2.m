X0=[0;0];

%Reaction Rates
%C(1)=\alpha_1 C(2)=\alpha_2 C(3)=\beta C(4)=\gamma
%C=[50;16;2.5;1];
C=[16;14;1;1];
%C=[];
% alpha_1=5; these parameters work
% alpha_2=5;
% beta=1;
% gamma=1;
alpha_1=C(1);
alpha_2=C(2);
beta=C(3);
gamma=C(4);
%Stoichiometry Matrix
S=[1 -1 0 0; 0 0 1 -1];
S_bis=[0 1 0 0; 0 0 0 1];



N_state=200;
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
W_star{1}=@(X,Y)(alpha_1/(1+(Y).^(beta)));
W_star{2}=@(X,Y)X(1);
W_star{3}=@(X,Y)(alpha_2/(1+(X).^(gamma)));
W_star{4}=@(X,Y)Y;
resample=1;
tf=5;
%%
[t,t_ob,delta,match,X_tot_prev,X,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);

%% FFSP1
[T,F,jump_times,E,Var,SD,Err_jump,E_tot,SD_tot,rho,eps] = FFSP_2_error(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);
%% FFSP2
[T_1,P,jump_times_1,E_1,Var_1,SD_1,Err_jump_1,E_tot_1,SD_tot_1,Err_tot_1] = FFSP_2_error_new(t_ob,Y,p0,C,State_Space,delta,S,S_bis,n1,W_star,MAK);



X_1_tot=X_tot_prev(n1,:); %exact trajectory of the hidden process 




%% 

f=figure;
f.Units='points';
f.OuterPosition=[10 10 1000 450];
subplot(1,2,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlabel('t (s)')
ylabel('Y(t)')
xlim([0 t_ob(end)])
title('Observation Process')
set(gca,'FontSize',20)  


freq1=30;
subplot(1,2,2)
stairs(t,X_1_tot(1,:),'b','LineWidth',2)
hold on
errorbar(t_ob',E,SD,'r--','LineWidth',3)
hold on
errorbar(t_ob,E_1,SD_1,'--','Color',[0.4660 0.6740 0.1880],'LineWidth',3.5)
hold off
title('Hidden Process')
xlim([0 t_ob(end)])
%ylim([0 15])
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
csvwrite('/Users/edambrosio/Desktop/PaperFigures/Toggle_Switch_Filter_1_Filter_2_bad_error',Tab);
