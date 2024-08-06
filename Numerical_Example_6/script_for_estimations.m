addpath(fullfile(pwd, 'Code_for_the_files'));


load('data45.mat')

mean_cell_intensity_open_loop

%% 
%genetic switch network
X0_1=[0;1;0];
C_1=[0.2,0.1,23,0.5];
W_star={};
MAK=1;

fluo=73;
t_1=timeAxis; %time of the single cell data analysis
%protein subnetwork 
S=[1 -1];
S_bis=[0 1];
c=[1,0.2]; %protein translation and degradation reaction rates 



n1=3;
n2=1;

N_state=200;

A=zeros(N_state+1,n1);
A(:,n1)=0:N_state;
A(:,1)=1;
B=ones(N_state+1,n1);
B(:,1)=0;
B(:,n1)=0:N_state;
State_Space=vertcat(A,B);


Shape_State_Space=size(State_Space);
State=Shape_State_Space(1);
[state,index]=intersect(State_Space,X0_1','rows');
p0=zeros(length(State_Space),1);
p0(index)=1;
S_2=[-1 1 0 0 0 0; 1 -1 0 0 0 0; 0 0 1 -1 0 0; 0 0 0 0 1 -1];
S_bis_2=[1 0 1 0 0 0; 0 1 0 0 0 0; 0 0 0 1 1 0; 0 0 0 0 0 1];


C_2=[C_1 c];

for m=1:30 %going through the cells data 
X_1=relevantCellTrajectory{1,m}.spotIntensity'./fluo;


Y0=0;
Y0bis=0;
Y=Y0;
t_ob=0;




for i=1:(length(t_1)-1)
    
    
    t0=t_1(i);
    tf=t_1(i+1);
    C=[c(1)*X_1(i),c(2)];
    
    [t,y] = Gillespie_General(Y0,S,S_bis,t0,tf,C,W_star,MAK);
    
    y=y(1:length(y)-1);
    t=t(1:length(t)-1);
    
    if length(y) > 1
        Y=[Y y(2:end)];
        t_ob=[t_ob t(2:end)];
        Y0=y(end);
    end
    
    
end







delta = diff(Y);
%[t_ob,Y,T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C_2,State_Space,S_2,S_bis_2,n1,W_star,MAK);
[T,F,jump_times,E_FSP,Var_FSP,SD_FSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0,C_2,State_Space,delta,S_2,S_bis_2,n1,W_star,MAK);

data_to_save{m}.t_ob=t_ob;
data_to_save{m}.Y=Y;
data_to_save{m}.E_FSP=E_FSP;
data_to_save{m}.Var_FSP=Var_FSP;
data_to_save{m}.SD_FSP=SD_FSP;
end

save path_of_cells_1 data_to_save;

%%
f=figure;
f.Units='points';
f.OuterPosition=[10 10 1000 950];
freq=1;
subplot(2,1,1)
stairs(t_ob,Y(1,:),'m','LineWidth',2)
xlim([0 t_ob(end)])
%xlabel("t (s) "+newline+"   ")
xlabel("t (min)")
ylabel('Y(t)')
title('Observation Process')
set(gca,'FontSize',20)  



subplot(2,1,2)
errorbar(t_ob([1 [1:freq:end]]),E_FSP([1 [1:freq:end]],3)',SD_FSP([1 [1:freq:end]],3),'r--','LineWidth',3)
hold on
stairs(t_1,X_1,'b','LineWidth',2)
hold off
title('Hidden Process (mRNA)')
%xlabel("t (s) "+newline+"   ")
xlabel(" t (min) ")
ylabel('Molecular Counts')
xlim([0 t_ob(end)])
legend('Exact Trajectory', 'FFSP','Location','northwest')
set(gca,'FontSize',20)  
