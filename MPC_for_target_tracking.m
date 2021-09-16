% This code is for the implementation of the MPC for the target - 
% tracking (without gimbal and without obstacles scenario) problem using CasADi v3.5.5.
clear all; clc; close all;

addpath('C:\Users\Dev Soni\Desktop\casadi')
import casadi.*

%% Defining the system ------------------------------------------------------------------------------------------

T = 0.2;
N = 15;

% Constrains of UAV

v_u_min = 14; v_u_max = 30;
omega_2_u_min = -pi/30; omega_2_u_max = pi/30;
omega_3_u_min = -pi/21; omega_3_u_max = pi/21;
theta_u_min = -0.2618; theta_u_max = 0.2618;
z_u_min = 75; z_u_max = 150;

% States without gimbal camera

x_u = SX.sym('x_u'); y_u = SX.sym('y_u'); z_u = SX.sym('z_u'); theta_u = SX.sym('theta_u'); psi_u = SX.sym('psi_u');   %states of the UAV
%x_t = SX.sym('x_t'); y_t = SX.sym('y_t'); theta_t = SX.sym('theta_t');                                               %states of the target

states_u = [x_u; y_u; z_u; theta_u; psi_u];      n_states_u = length(states_u); %UAV
%states_t = [x_t; y_t; theta_t];    n_states_t = length(states_t); %target

% Controls

v_u = SX.sym('v_u'); omega_2_u = SX.sym('omega_2_u'); omega_3_u = SX.sym('omega_3_u');        % UAV cotrols parameters
%v_t = SX.sym('v_t'); omega_t = SX.sym('omega_3');                                             % control Parameters of target

controls = [v_u; omega_2_u; omega_3_u]; n_controls = length(controls);
rhs = [v_u*cos(psi_u)*cos(theta_u); v_u*sin(psi_u)*cos(theta_u); v_u*sin(theta_u);...
    omega_2_u; omega_3_u];                     %systems r.h.s

f = Function('f',{states_u,controls}, {rhs});                                                 % Nonlinear Mapping Function f(x,u)
U = SX.sym('U',n_controls, N);                                                                % Desition Variables 
P = SX.sym('P',n_states_u + 2);                                                        
% This consists of initial states od UAV 1-5 and reference states 6-7(reference states are target's states)

X = SX.sym('X',n_states_u,(N+1));
% Consists of future predicted states of UAV

%% Filling the defined sysem parameters --------------------------------------------------------------------------

X(:,1) = P(1:5) % initial state
for k = 1:N
    st = X(:,k); con = U(:,k);
    f_value = f(st,con);
    st_next = st + T*f_value;
    X(:,k+1) = st_next;
end

ff = Function('ff',{U,P},{X});

%% Objective function --------------------------------------------------------------------------------------------

obj = [];  % object function
g = [];    % constrain vector

w1 = 0.1; %w2 = 1; % weighting coefficents each for one function w1 -- f1 and w2 -- f2

for k = 1:N
    x = states_u(1:1); y = states_u(2:1);
    %obj = obj + sqrt((x - target's x)^2 + (y - target's y)^2);
    
end


%% simulation of the target and UAV's distance ---------------------------------------------------------------------

t = 0:T:10 %TIme Span

% Initial Conditions
x_u_0 = 90;
y_u_0 = 150;     
z_u_0 = 80;     % for UAV

states_u(:,1) = [x_u_0; y_u_0;z_u_0;0;0];

x_t_0 = 100;
y_t_0 = 150;
theta_t_0 = 0; % for target

eta0 = [x_t_0; y_t_0; theta_t_0];

states_t(:,1) = eta0;

% Loop started from here

for i=1:length(t)
    
    shi = states_t(3,i); % Orientation of target
    v_t = 12;  % linear velocity of target
    omega_t = 1; % Angular velocity of target
    z(1,i)=0; %target is on ground
    
    eta_dot(:,i) = [v_t*cos(shi); v_t*sin(shi); omega_t];
    
    states_t(:,i+1) = states_t(:,i) + T * eta_dot(:,i);
    
    cost_value(1,i) = sqrt((x_u_0-states_t(1,i))^2 + (y_u_0-states_t(2,i))^2);
end


subplot(2,3,5) %plotting the cost values
plot(t,cost_value(1,1:i));
xlabel('T[s]');
ylabel('cost value');
grid on;
    
    
    
