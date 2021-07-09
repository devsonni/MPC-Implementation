% This code is for the implementation of the MPC for the target - 
% tracking problem using CasADi v3.5.5.

addpath('C:\Users\Dev Soni\Desktop\casadi')
import casadi.*

%% Defining the system 

T = 0.2;
N = 15;

% Constrains of UAV

v_u_min = 14; v_u_max = 30;
omega_2_u_min = -pi/30; omega_2_u_max = pi/30;
omega_3_u_min = -pi/21; omega_3_u_max = pi/21;
theta_u_min = -0.2618; theta_u_max = 0.2618;
z_u_min = 75; z_u_max = 150;

% States without gimbal camera

x_u = SX.sym('x_u'); y_u = SX.sym('y_u'); z_u = SX.sym('z_u'); theta_u = SX.sym('theta_u'); psi_u = SX.sym('psi_u'); %states of the UAV
x_t = SX.sym('x_t'); y_t = SX.sym('y_t'); theta_t = SX.sym('theta_t');                                               %states of the target

states = [x_u; y_u; z_u; theta_u; psi_u; x_t; y_t; theta_t];      n_states = length(states);

% Controls

v_u = SX.sym('v_u'); omega_2_u = SX.sym('omega_2_u'); omega_3_u = SX.sym('omega_3_u');      % UAV cotrols
v_t = SX.sym('v_t'); omega_t = SX.sym('omega_3');                                           % Parameters of target
controls = [v_u; omega_2_u; omega_3_u]; n_controls = length(controls);
rhs = [v_u*cos(psi_u)*cos(theta_u); v_u*sin(psi_u)*cos(theta_u); v_u*sin(theta_u);...
    omega_2_u; omega_3_u; v_t*cos(theta_t); v_t*sin(theta_t); omega_t];                     %systems r.h.s

f = Function('f',{states,controls}, {rhs});                                                 % Nonlinear Mapping Function f(x,u)
U = SX.sym('U',n_controls, N);                                                              % Desition Variables 
P = SX.sym('P',n_states + n_states);                                                       
% This consists of initial states 0-7 and reference states 8-15

X = SX.sym('X',n_states,(N+1));
% Consists of future predicted states

%% Filling the defined sysem parameters

X(:,1) = P(1:8) % initial state
for k = 1:N
    st = X(:,k); con = U(:,k);
    f_value = f(st,con);
    st_next = st + T*f_value;
    X(:,k+1) = st_next;
end

ff = Function('ff',{U,P},{X});
    
    
    
    
    
    
    
