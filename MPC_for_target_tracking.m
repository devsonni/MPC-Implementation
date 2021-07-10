% This code is for the implementation of the MPC for the target - 
% tracking (without gimbal and without obstacles scenario) problem using CasADi v3.5.5.

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

%% Filling the defined sysem parameters --------------------------------------------------------------------------

X(:,1) = P(1:8) % initial state
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

w1 = 0.1; w2 = 1; % weighting coefficents each for one function w1 -- f1 and w2 -- f2

for k=1:N
    
    x_u = states(1:1);   y_u = states(2:1); x_t = states(6:1); y_t = states(7:1);
    a = z_u*tan(VFOV/2); b = z_u*tan(HFOV/2);                                          % defining field of VFOV and HFOV is rest to done
    A_k = (cos(psi)/a)^2 + (sin(psi)/b)^2 ; B_k = 2*cos(psi)*sin(psi)*(1/a^2 - 1/b^2) ; C_k = (sin(psi)/a)^2 + (cos(psi)/b)^2;  
    obj = w1*sqrt((x_u - x_t)^2+(y_u - y_t)^2) + w2*((A_k*(x_u - x_t)^2 + B_k*(y_t - y_u)*(x_t - x_u) + C_k*(y_t - y_u)^2)-1); % J(X,U) objective function
    
end

    
    
    
    
    
    
    
