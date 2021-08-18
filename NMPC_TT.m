% This code is for the implementation of the MPC for the target - 
% tracking (without gimbal and without obstacles scenario) problem using CasADi v3.5.5.
clear all; clc; close all;

addpath('C:\Users\Dev Soni\Desktop\casadi')
import casadi.*

%% Defining the system of UAV ------------------------------------------------------------------------------------------

T = 0.2;
N = 15;

% Constrains of UAV

v_u_min = 14; v_u_max = 30;
omega_2_u_min = -pi/30; omega_2_u_max = pi/30;
omega_3_u_min = -pi/21; omega_3_u_max = pi/21;
theta_u_min = -0.2618; theta_u_max = 0.2618;
z_u_min = 75; z_u_max = 150;

% States of UAV with-out gimbal camera

x_u = SX.sym('x_u'); y_u = SX.sym('y_u'); z_u = SX.sym('z_u'); theta_u = SX.sym('theta_u'); psi_u = SX.sym('psi_u');   %states of the UAV
states_u = [x_u; y_u; z_u; theta_u; psi_u];      n_states_u = length(states_u);               %UAV

% Controls of UAV that will find by NMPC

v_u = SX.sym('v_u'); omega_2_u = SX.sym('omega_2_u'); omega_3_u = SX.sym('omega_3_u');        % UAV cotrols parameters

controls_u = [v_u; omega_2_u; omega_3_u]; n_controls = length(controls_u);
rhs_u = [v_u*cos(psi_u)*cos(theta_u); v_u*sin(psi_u)*cos(theta_u); v_u*sin(theta_u);...
    omega_2_u; omega_3_u];                                                                    %systems r.h.s

f_u = Function('f_u',{states_u,controls_u}, {rhs_u});                                         % Nonlinear Mapping Function f(x,u)
U = SX.sym('U',n_controls, N);                                                                % Desition Variables 
P = SX.sym('P',n_states_u + 2);                                                        
% This consists of initial states od UAV 1-5 and reference states 6-7(reference states are target's states)

X = SX.sym('X',n_states_u,(N+1));
% Consists of future predicted states of UAV

% Filling the defined sysem parameters of UAV --------------------------------------------------------------------------

X(:,1) = P(1:5) % initial state
for k = 1:N
    st = X(:,k); con = U(:,k);
    f_value = f_u(st,con);
    st_next = st + T*f_value;
    X(:,k+1) = st_next;
end

ff = Function('ff',{U,P},{X});

%% Defining the system of Target ------------------------------------------------------------------------------------------

%states of target
x_t = SX.sym('x_t'); y_t = SX.sym('y_t'); theta_t = SX.sym('theta_t');                                               
states_t = [x_t; y_t; theta_t];    n_states_t = length(states_t); 

%Controls of target
v_t = SX.sym('v_t'); omega_t = SX.sym('omega_t');  

controls_t = [v_t; omega_t];
rhs_t = [v_t*cos(theta_t); v_t*sin(theta_t); omega_t];

f_t = Function('f_t',{states_t,controls_t}, {rhs_t}); 

Y = SX.sym('Y',n_states_t,N+1);

x_t_0 = 100;
y_t_0 = 150;
theta_t_0 = 0; 
v_t = 12;
omega_t = 1;

Y(:,1) = [x_t_0; y_t_0; theta_t_0];

for k=1:N
    st = Y(:,k); con = [v_t; omega_t];
    f_value = f_t(st,con);
    st_next = st + T*f_value;
    Y(:,k+1) = st_next;   
    reference = Y(1:2,1:N);
end

obj = 0; %objective function
g = [];  %constrains

for k=1:N
    stt = X(1:2,1:N); reference(1:2,k);
    obj = obj + sqrt((stt(1,k)-reference(1,k))^2 + (stt(2,k)-reference(2,k))^2);
end

 
