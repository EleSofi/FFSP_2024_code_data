%reaction rate parameters

Kr=100;

gamma_r=10;
kp=5;
%gamma_p=2;
%gamma_fb=6;
omega=100;
kr=Kr/omega;
N=1;
n1=1;

MAK=1;
X0=[0;0];

W_star=[];
C=[Kr,gamma_r,kp];
S=[1 -1 0; 0 0 1];
S_bis=[0 1 1; 0 0 0];
resample=1;
N_tot=10000;
tf=1;

%M_tot=x(1,:)+(sqrt(1/omega)*m);
    

N_state=500;

Hidden_Species=Hidden_State(0: N_state);
Shape_State_Space_Hidden=size(Hidden_Species);
State_Hidden=Shape_State_Space_Hidden(1);
[c,index]=intersect(Hidden_Species,X0(1:n1)','rows');
p0_1=zeros(State_Hidden,1);
p0_1(index)=1;

%% gillespie simulation for generating the hidden X_tot_prev process (time of the whole network t_1) and observation process Y 
[t_1,t_ob,delta_1,match,X_tot_prev,X_1,Y] = Gillespie_General_MAK(X0,S,S_bis,tf,C,n1,W_star,MAK);
Y_omega=Y./omega;
%% particle filter 
[V_tot_1,w_tot,V_jump,w_jump,match_1,match_2,resampling,pt,E_PF,Var_PF,SD_PF] = particle_filter_1(t_ob,Y,p0_1,C,tf,N_tot,Hidden_Species,S,S_bis,resample,n1,W_star,MAK);

%%
%scaled process 
X_1=X_1./omega;
X_1_tot=X_tot_prev(n1,:)./omega;
%time flow of the process

%% FFSP filter 
[T,F,jump_times,E_FFSP,Var_FFSP,SD_FFSP,Err_jump,E_tot,SD_tot,rho] = FFSP_2(t_ob,Y,p0_1,C,Hidden_Species,delta_1,S,S_bis,n1,W_star,MAK);
% E_FSP=E_FSP./omega;
% SD_FSP=SD_FSP./omega;
% X_1=X_1./omega; %exact trajectory of the hidden process at the observation process jump times
% X_1_tot=X_tot_prev(n1,:)./omega; %exact trajectory of the hidden process 


    
%%
%Kalman filter implementation on top of LNA


% Initialization
N = length(t_ob); % Number of time points
x = zeros(2, N);
x(1, 1) = 10^-3;  % Initial condition for x1
x(2, 1) = 0;      % Initial condition for x2, adjust if different initial condition is needed
mu = zeros(1, N);
sigma = zeros(1, N);

mu(1) = 0;        % Initial condition for mu
sigma(1) = 0;     % Initial condition for sigma

% Constants initialization, assuming all are defined elsewhere in your code
for i = 1:N-1
    delta = t_ob(i+1) - t_ob(i);
    
    % Update system states
    x(1, i+1) = x(1, i) + delta * (kr - gamma_r * x(1, i));
    x(2, i+1) = x(2, i) + delta * (kp * x(1, i));
    
    % Conditional mean evolution equation
    increment_Y = Y_omega(i+1) - Y_omega(i);
    increment_x2 = x(2, i+1) - x(2, i);
    mu(i+1) = mu(i) + delta * (-gamma_r * mu(i)) + ...
              (sigma(i) * x(1, i)^-1) * (sqrt(omega) * (increment_Y - increment_x2) - kp * mu(i) * delta);
    
    % Conditional variance evolution equation (Kushner equation)
    D = [sqrt(kr); -sqrt(gamma_r * x(1, i))];
    sigma(i+1) = sigma(i) + delta * (2 * (-gamma_r) * sigma(i) + ...
                    D' * D - (sigma(i)^2) * (kp / x(1, i)));
end

    
   


E_LNA=x(1,:)+(mu/sqrt(omega));
Sigma=sigma./omega;
SD_LNA=sqrt(Sigma);



%% plotting the filters results with the observation trajectory 

f = figure;
f.Units = 'points';
f.OuterPosition = [10 10 1000 450];
freq = 2;

% Observation Process Plot
subplot(1, 2, 1)
stairs(t_ob, Y, 'm', 'LineWidth', 3)
xlabel('t (s)')
ylabel('Y(t)')
title('Observation Process')
xlim([0 t_ob(end)])
set(gca, 'FontSize', 20)
ax1 = gca;
ax1.XAxis.Exponent = 0;
ax1.YAxis.Exponent = 2;

% Hidden Process Plot
subplot(1, 2, 2)
stairs(t_1, X_tot_prev(1, :), 'b', 'LineWidth', 3)
hold on
errorbar(t_ob([1 [1:freq:end]]), E_PF(([1 [1:freq:end]]), 1)', SD_PF(([1 [1:freq:end]]), 1), '-.', 'Color', "#EDB120", 'LineWidth', 1.5)
errorbar(t_ob, E_LNA * omega, SD_LNA * omega, '--', 'LineWidth', 1.5, 'Color','c')
errorbar(t_ob([1 [1:freq:end]]), E_FFSP([1 [1:freq:end]]), SD_FFSP([1 [1:freq:end]]), 'r--', 'LineWidth', 1.5)
xlim([0 t_ob(end)])
xlabel('t (s)')
ylabel('Molecular Counts')
legend('Exact', 'BPF', 'Kalman', 'FFSP', 'Location', 'northwest')
title('Hidden Process')
set(gca, 'FontSize', 20)
ax2 = gca;
ax2.XAxis.Exponent = 0;
ax2.YAxis.Exponent = 2;
