clear all; clc; close all;
addpath('C:\Users\Dev Soni\Desktop\casadi')
import casadi.*
addpath('C:\Users\Dev Soni\Desktop\npy-matlab-master\npy-matlab')
import npy-matlab.*

x0 = [99;150;80;0;0;0;0;0]; %initial location of UAV and gimbal
xs = [100;150;0]; %reference for the UAV which is target
%Parameters of Obstacle
x_o_1 = 500;
y_o_1 = 20;
x_o_2 = 1600;
y_o_2 = 197;
x_o_3 = 130;
y_o_3 = 670;
obs_r = 100;
UAV_r = 5;

arr = unzip('500epi_095alpha_085gamma_03min_0074decay.npz');
max_step_size = readNPY('max_step_size.npy');
max_steps = max_step_size;
qtable = readNPY('qtable.npy');
total_episodes = readNPY('total_episodes.npy');

% Parameters for plotting
rewardarr = zeros(total_episodes, max_steps);
errorarr = zeros(total_episodes, max_steps);
erroravg = zeros(total_episodes);
rewardavg = zeros(total_episodes);
maxreward = zeros(max_steps);
minerrorarr = zeros(max_steps);
action_str = zeros(max_steps + 1, 2);
UAV_RL = zeros(8,max_step_size+1);
UAV_W_RL = zeros(8,max_step_size+1);
targetarr = zeros(3,max_step_size+1);
error_w_rl = zeros(max_steps);
FOV_C_RL = zeros(max_step_size, 2);
FOV_C_W_RL = zeros(max_step_size, 2);

%Filling optimal policy stuff
sc = 0;
error = 0;
UAV_RL(:, 1) = x0;
targetarr(:,1)= xs;
for i=1:max_step_size
    index = MAXINDEX(i,qtable);
    [Error, x0, sc, xs, a_p, b_p, FOV_X, FOV_Y] = MPC(index(1)-1, index(2)-1, x0, xs, sc);
    sc
    FOV_C_RL(i, 1) = FOV_X;
    FOV_C_RL(i, 2) = FOV_Y;
    UAV_RL(:, i+1) = x0;
    targetarr(:,i+1)= xs;
    maxreward(i) = 1/Error;
    minerrorarr(i) = Error;
    error = error + Error;
end
[x_e_rl, y_e_rl] = ellipse(a_p, b_p, FOV_X, FOV_Y);
Error_RL = error

%Resetting Initial location of UAV and target
x0 = [99;150;80;0;0;0;0;0]; %initial location of UAV and gimbal
xs = [100;150;0]; %reference for the UAV which is target

%filling without RL stuff
sc = 0;
error = 0;
UAV_W_RL(:, 1) = x0;
for i=1:max_step_size
    [Error, x0, sc, xs, a_p, b_p, FOV_X, FOV_Y] = MPC(1, 2, x0, xs, sc);
    sc
    FOV_C_W_RL(i, 1) = FOV_X;
    FOV_C_W_RL(i, 2) = FOV_Y;
    UAV_W_RL(:, i+1) = x0;
    error_w_rl(i) = Error;
    error = error + Error;
end
[x_e_w_rl, y_e_w_rl] = ellipse(a_p, b_p, FOV_X, FOV_Y);
Error_w_RL = error


%Plotting starts from here

% just for  plotting 
ss1 = [];
ss1(1,1:(max_steps+100)) = 0;
[X,Y,Z] = cylinder(obs_r);
h = 150;
Z = Z*h;

ss2 = zeros(max_steps);
for i=1:max_steps
    ss2(i) = i;
end


figure
plot3(UAV_RL(1,1:max_steps),UAV_RL(2,1:max_steps),UAV_RL(3,1:max_steps),'b-',...
    UAV_W_RL(1,1:max_steps),UAV_W_RL(2,1:max_steps),UAV_W_RL(3,1:max_steps),'r-',...
    targetarr(1,1:max_steps),targetarr(2,1:max_steps),ss1(1,1:max_steps), 'r--',...
    x_e_rl, y_e_rl, ss1(1,1:101), 'g-',...
    x_e_w_rl, y_e_w_rl, ss1(1,1:101),'y-', 'LineWidth', 3);
hold on;
surf(X+x_o_1,Y+y_o_1,Z);
hold on;
surf(X+x_o_2,Y+y_o_2,Z);
hold on;
surf(X+x_o_3,Y+y_o_3,Z);
xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
legend('UAV With RL', 'UAV Without RL', 'Target', 'FOV With RL', 'FOV Without RL', 'Obstacles');
xlim([-300 2100]);
ylim([-400 1100]);
zlim([0 500]);
grid on;

figure
subplot(1,2,1)
plot(ss2, minerrorarr, 'b-', ss2, error_w_rl, 'r-', 'LineWidth', 2);
xlabel('Steps')
ylabel('Error [m]');
legend('Error With RL', 'Error Without RL');
title('Error Trajectories');

subplot(1,2,2)
x = {'Error Of FOV With RL', 'Error Of FOV Without RL'};
y = [Error_RL, Error_w_RL];
barh(y)
set(gca,'yticklabel', x)
title('Sum Of FOV-E With And Without RL');

