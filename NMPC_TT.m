% This code is for the implementation of the MPC for the target - 
% tracking (without gimbal and without obstacles scenario) problem using CasADi v3.5.5.
clear all; clc; close all;

addpath('C:\Users\Dev Soni\Desktop\casadi')
import casadi.*

%% Defining the system of UAV ------------------------------------------------------------------------------------------

T = 0.2;
N = 15;

% Constrains of UAV

  %input constrains
v_u_min = 14; v_u_max = 30;
omega_2_u_min = -pi/30; omega_2_u_max = pi/30;
omega_3_u_min = -pi/21; omega_3_u_max = pi/21;

  %states constrains
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
P = SX.sym('P',n_states_u + 3);                                                        
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

%{
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
%}

obj = 0; %objective function
g = [];  %constrains of pitch angle theta

%{
for k=1:N
    stt = X(1:2,1:N); reference(1:2,k);
    obj = obj + sqrt((stt(1,k)-reference(1,k))^2 + (stt(2,k)-reference(2,k))^2);
end
%}

for k=1:N
    stt = X(1:2,1:N); 
    obj = obj + sqrt((stt(1,k)-P(6))^2 + (stt(2,k)-P(7))^2);
end

%compute the constrains
 for k=1:N+1
     g = [g, X(3,k)]; %limit on height
     g = [g, X(4,k)]; %limit on pitch angle theta
 end
 
% make the decision variables one column vector
OPT_variables = reshape(U,3*N,1);
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
args.lbg(1,1:2:32) = 75;     args.lbg(1,2:2:32) = -0.2618;    % lower bound of height and theta 
args.ubg(1,1:2:32) = 150;    args.ubg(1,2:2:32) =  0.2618;    % upper bound of height and theta 

% input constraints
args.lbx(1:3:3*N-1,1) = v_u_min; args.lbx(2:3:3*N,1) = omega_2_u_min; args.lbx(3:3:3*N,1) = omega_3_u_min;
args.ubx(1:3:3*N-1,1) = v_u_max; args.ubx(2:3:3*N,1) = omega_2_u_max; args.ubx(3:3:3*N,1) = omega_3_u_max;

%% Simulation starts from here----------------------------------------------------------------------------------------------

t0 = 0;
x0 = [90;150;80;0;0]; %initial location of UAV
xs = [100;150;0]; %reference for the UAV which is target
xx(:,1) = x0; %storing history of location the UAV
t(1) = t0;
u0 = zeros(N,3); %initial control of UAV
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
    args.x0 = reshape(u0',3*N,1); % initial value of the optimization variables
    %tic
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
            'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);    
    %toc
    u = reshape(full(sol.x)',3,N)';
    ff_value = ff(u',args.p); % compute OPTIMAL solution TRAJECTORY
    xx1(:,1:5,mpciter+1)= full(ff_value)';
    xss(:,1:3,mpciter+1) = full(xs)'; 
    
    u_cl= [u_cl ; u(1,:)];
    t(mpciter+1) = t0;
    [t0, x0, u0, xs] = shift1(T, t0, x0, u, f_u, f_t, xs); % get the initialization of the next optimization step
    
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

%{ 
subplot(1,2,1)
plot3(ss(1,1:100),ss(2,1:100), ss1(1,1:100));
xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
grid on;
%}
figure
plot3(ss(1,1:100),ss(2,1:100), ss1(1,1:100),xx(1,1:100),xx(2,1:100), xx(3,1:100));
xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
legend('Target','UAV');
grid on;


%% Simulation --------------------------------------------------------------------------------------------------------------------------
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

myWriter = VideoWriter('Tracking4');
myWriter.FrameRate = 20;

open(myWriter);
writeVideo(myWriter, movieVector);
close(myWriter);
