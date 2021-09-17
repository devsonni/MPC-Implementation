% This code is for the implementation of the MPC for the target - 
% tracking (without gimbal and without obstacles scenario) problem using CasADi v3.5.5.
clear all; clc; close all;

addpath('C:\Users\Dev Soni\Desktop\casadi')
import casadi.*

%% Defining the system of UAV ------------------------------------------------------------------------------------------

T = 0.2;
N = 15;

% Constrains of UAV with gimbal

  % input constrains of UAV
v_u_min = 14; v_u_max = 30;
omega_2_u_min = -pi/30; omega_2_u_max = pi/30;
omega_3_u_min = -pi/21; omega_3_u_max = pi/21;

  % input constrains of gimbal
omega_1_g_min = -pi/30; omega_1_g_max = pi/30;
omega_2_g_min = -pi/30; omega_2_g_max = pi/30;
omega_3_g_min = -pi/30; omega_3_g_max = pi/30;
  
  % states constrains of UAV
theta_u_min = -0.2618; theta_u_max = 0.2618;
z_u_min = 75; z_u_max = 150;

  % states constrains of gimbal
phi_g_min = -pi/6; phi_g_max = pi/6;
theta_g_min = -pi/6; theta_g_max = pi/6;
shi_g_min = -pi/2; shi_g_max = pi/2;

% States of UAV with gimbal camera

x_u = SX.sym('x_u'); y_u = SX.sym('y_u'); z_u = SX.sym('z_u'); theta_u = SX.sym('theta_u'); psi_u = SX.sym('psi_u');   %states of the UAV
phi_g = SX.sym('phi_g'); shi_g = SX.sym('shi_g'); theta_g = SX.sym('theta_g');                                         %states of gimbal

states_u = [x_u; y_u; z_u; theta_u; psi_u; phi_g; theta_g; shi_g];      n_states_u = length(states_u);               %UAV with gimbal

% Controls of UAV that will find by NMPC

v_u = SX.sym('v_u'); omega_2_u = SX.sym('omega_2_u'); omega_3_u = SX.sym('omega_3_u');              % UAV cotrols parameters
omega_1_g = SX.sym('omega_1_g'); omega_2_g = SX.sym('omega_2_g'); omega_3_g = SX.sym('omega_3_g');  % Gimbal control parameters

controls_u = [v_u; omega_2_u; omega_3_u; omega_1_g; omega_2_g; omega_3_g]; n_controls = length(controls_u);
rhs_u = [v_u*cos(psi_u)*cos(theta_u); v_u*sin(psi_u)*cos(theta_u); v_u*sin(theta_u);...
    omega_2_u; omega_3_u; omega_1_g; omega_2_g; omega_3_g];                                         %systems r.h.s

f_u = Function('f_u',{states_u,controls_u}, {rhs_u});                                         % Nonlinear Mapping Function f(x,u)
U = SX.sym('U',n_controls, N);                                                                % Desition Variables 
P = SX.sym('P',n_states_u + 3);                                                        
% This consists of initial states of UAV with gimbal 1-8 and reference states 9-11 (reference states are target's states)

X = SX.sym('X',n_states_u,(N+1));
% Consists of future predicted states of UAV

% Filling the defined sysem parameters of UAV --------------------------------------------------------------------------

X(:,1) = P(1:8) % initial state
for k = 1:N
    st = X(:,k); con = U(:,k);
    f_value = f_u(st,con);
    st_next = st + T*f_value;
    X(:,k+1) = st_next;
end

ff = Function('ff',{U,P},{X});

%% Objective function-----------------------------------------------------------------------------------------------------

obj = 0; %objective function
g = [];  %constrains of pitch angle theta

for k=1:N
    stt = X(1:8,1:N); A= SX.sym('A',N); B= SX.sym('B',N); C= SX.sym('C',N); X_E= SX.sym('X_E',N); Y_E= SX.sym('Y_E',N);
    VFOV = 5;
    HFOV = 10;
    w1 = 0.1;
    w2 = 1;
    a = (stt(3,k)*(tan(stt(7,k)+VFOV/2)) - stt(3,k)*tan(stt(7,k)-VFOV/2))/2;
    b = (stt(3,k)*(tan(stt(6,k)+HFOV/2)) - stt(3,k)*tan(stt(6,k)-HFOV/2))/2;
    A(k) = ((cos(stt(8,k)))^2)/a^2 + ((sin(stt(8,k)))^2)/b^2;
    B(k) = 2*cos(stt(8,k))*sin(stt(8,k))*((1/a^2)-(1/b^2));
    C(k) = ((sin(stt(8,k)))^2)/a^2 + ((cos(stt(8,k)))^2)/b^2;
    X_E(k) = stt(1,k) + a + stt(3,k)*(tan(stt(7,k)-VFOV/2));
    Y_E(k) = stt(2,k) + b + stt(3,k)*(tan(stt(6,k)-HFOV/2));
    obj = obj + w1*sqrt((stt(1,k)-P(9))^2 + (stt(2,k)-P(10))^2) + w2*((A(k)*(P(9)-X_E(k))^2 +...
        B(k)*(P(10)-Y_E(k))*(P(9)-X_E(k)) + C(k)*(P(10)-Y_E(k))^2)-1);
end

%compute the constrains
 for k=1:N+1
     g = [g, X(3,k)]; %limit on height
     g = [g, X(4,k)]; %limit on pitch angle theta
     g = [g, X(6,k)]; %limit on gimbal angle phi
     g = [g, X(7,k)]; %limit on gimbal angle theta
     g = [g, X(8,k)]; %limit on gimbal angle shi
 end
 
% make the decision variables one column vector
OPT_variables = reshape(U,6*N,1);
nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', P);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;

% inequality constraints (state constraints)
args.lbg(1,1:5:80) = 75;          args.ubg(1,1:5:80) = 150;             % bounds on height  
args.lbg(1,2:5:80) = -0.2618;     args.ubg(1,2:5:80) =  0.2618;         % bounds on theta of UAV 
args.lbg(1,3:5:80) = phi_g_min;   args.ubg(1,3:5:80) = phi_g_max;       % bounds on gimbal angle phi
args.lbg(1,4:5:80) = theta_g_min; args.ubg(1,4:5:80) = theta_g_max;     % bounds on gimbal angle theta
args.lbg(1,5:5:80) = shi_g_min;   args.ubg(1,5:5:80) = shi_g_max;       % bounds on gimbal angle shi

% input constraints
args.lbx(1:6:6*N-1,1) = v_u_min;  args.lbx(2:6:6*N,1) = omega_2_u_min;  args.lbx(3:6:6*N,1) = omega_3_u_min;  % UAV's input bounds
args.ubx(1:6:6*N-1,1) = v_u_max;  args.ubx(2:6:6*N,1) = omega_2_u_max;  args.ubx(3:6:6*N,1) = omega_3_u_max;

args.lbx(4:6:6*N,1) = omega_1_g_min;     args.ubx(4:6:6*N,1) = omega_1_g_max;      % gimbal's input bounds
args.lbx(5:6:6*N,1) = omega_2_g_min;     args.ubx(5:6:6*N,1) = omega_2_g_max; 
args.lbx(6:6:6*N,1) = omega_3_g_min;     args.ubx(6:6:6*N,1) = omega_3_g_max;

%% Simulation starts from here----------------------------------------------------------------------------------------------

t0 = 0;
x0 = [90;150;80;0;0;0;0;0]; %initial location of UAV and gimbal
xs = [100;150;0]; %reference for the UAV which is target
xx(:,1) = x0; %storing history of location the UAV
t(1) = t0;
u0 = zeros(N,6); %initial control of UAV
sim_time = 20;

% NMPC starts form here

mpciter = 0;
xx1 = [];
u_cl = [];
xss = [];
ss = [];
ss(:,1) = xs;

main_loop = tic;
while (mpciter < sim_time/T)
    args.p = [x0;xs];
    args.x0 = reshape(u0',6*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',6,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:8,mpciter+1)= full(ff_value)';
    xss(:,1:3,mpciter+1) = full(xs)'; 
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0, xs] = shift1(T, t0, x0, u, f_u, xs); % get the initialization of the next optimization step
    
    ss(:,mpciter+2) = xs;
    xx(:,mpciter+2) = x0;  
    mpciter
    mpciter = mpciter + 1;

end
main_loop_time = toc(main_loop)

%% Plotting stuff----------------------------------------------------------------------------------------------------------------------

% just for  plotting 
ss1 = [];
ss1(1,1:100) = 0;


figure
plot3(ss(1,1:100),ss(2,1:100), ss1(1,1:100),xx(1,1:100),xx(2,1:100), xx(3,1:100));
xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
legend('Target','UAV');
grid on;
 

%{
% Code for saving video file of simulation

figh = figure
for i=1:100
    hold on;
    plot3(ss(1,i),ss(2,i), ss1(1,i),xx(1,i),xx(2,i), xx(3,i),'go', 'LineWidth', 4, 'MarkerSize', 1);
    plot3(ss(1,i),ss(2,i), ss1(1,i), 'bo', 'LineWidth', 4, 'MarkerSize', 1);
    xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
    legend('Target','UAV');
    view([30+i 35]);
    grid on;
    hold off;
    drawnow;
    
    movieVector(i) = getframe(figh, [10 10 520 400]);
    
end

myWriter = VideoWriter('Tracking6');
myWriter.FrameRate = 10;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
%}