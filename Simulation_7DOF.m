clear 
clc 
clear all 
%% Simulation 
% l = [1;1];  
% m = [1;1]; 
%q = [-pi/2;0]; 
%dq = [-pi/2;0]; 
%ddq = [-pi/2;0];
% i = [1;1]  %Trägheiten in einem Vektor zusammengefasst
% gp = 9.81; 
% dq = [1;0];
% dq1 = 1;
% dq2 = 1;
% 
% M = M_n2(q,l,m,i)
% C = C_n2(q,l,m,dq)
% g = g_n2(q,l,m,gp)
%% Start simulation in Simulink  
cxyz = zeros(1,21);
inertia_i = [0.1,0,0,0,0.09,0,0,0,0.02,0.05,0,0,0,0.018,0,0,0,0.044,0.08,0,0,0,0.075,0,0,0,0.01,0.03,0,0,0,0.01,0,0,0,0.029,0.02,0,0,0,0.018,0,0,0,0.05,0.005,0,0,0,0.0036,0,0,0,0.0047,0.001,0,0,0,0.001,0,0,0,0.001]'; 
%initial_condition_q7 = [0,pi/2,0,pi/2,0,-pi/2,0];
time = 50; %Duration of the simulation in seconds  
initial_condition_q7 = [0,0,0,0,0,0,0]; %Starting angles of the joints
simOut = sim("C:\Users\djama\OneDrive\Desktop\Djamal\Bachelorarbeit\Matlab\7DOF\Differntial_Equation7DOF.slx",[0 time]);

%% Measured Values from Simulink 
time = squeeze(0:0.1:time);
%time = linspace(0,10,101);
  gp = 9.81;



%ideal, without noise 
Tau_n7 = squeeze(simOut.Tau_sim7.data);
qSim7 = squeeze(simOut.qSim7.data);     
dqSim7 = squeeze(simOut.dqSim7.data); 
ddqSim7 = squeeze(simOut.ddqSim7.data);
Tau_sim7 = [];
% Y_sim = [];
% Y_sim7N = []; 

%build stacked regressor matrix and Tau vector
Tau_sim7 = reshape(Tau_n7, length(Tau_n7)*7,1);
% for i = 1:length(time)
%     Y_sim = [Y_sim;Y_n7(qSim7(:,i),dqSim7(:,i),ddqSim7(:,i),gp)];
% end 

%With noise in output 
Tau_sim7N = Tau_sim7 + rand(size(Tau_sim7)) * pi/1800; 

%With noise in outputs and inputs  
% qr = rand(2,length(time)) * pi/1800;    %White noise with amplitude = 0,1°, frequency = 1/2.75
% dqr = rand(2,length(time)) * pi/1800 * 2 * pi * 1/2.75 ; 
% ddqr = rand(2,length(time)) * pi/1800 * (2 * pi * 1/2.75)^2 ;
% 
% qSim7N = qSim7 + qr;    %Add noise to the ideal simulation 
% dqSim7N = dqSim7 + dqr; 
% ddqSim7N = ddqSim7 + ddqr; 

% for i = 1:length(time)
%     Y_sim7N = [Y_sim7N;Y_n7(qSim7N(:,i),dqSim7N(:,i),ddqSim7N(:,i),gp)];
% end 
% 
% 
% p1 = lsqr(Y_sim,Tau_sim); 
% p2 = lsqr(Y_sim,Tau_sim2); 
% p3 = lsqr(Y_sim1,Tau_sim2); 
% 
% c1 = cond(Y_sim'*Y_sim); 
% c2 = cond(Y_sim1'*Y_sim1); 


