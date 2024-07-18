% Target state pridiction model
clear all; clc; close all;

%Simulation Parameters
dt = 0.1; %tim-size
ts = 10; %Simulation time
t = 0:dt:ts; %Time span

%initial location of the target 
x0 = 100;
y0 = 150;
shi0 = 0;

%initial location of the UAV (90,150,80)
x_u0 = 90;
y_u0 = 150;
z_u0 = 80;

%intermediate matrix
eta0 = [x0; y0; shi0];
eta(:,1) = eta0;


for i=1:length(t)
    
    shi = eta(3,i); % Orientation of target
    v_t = 12;  % linear velocity of target
    omega_t = 1; % Angular velocity of target
    z(1,i)=0; %target on ground
    
    eta_dot(:,i) = [v_t*cos(shi); v_t*sin(shi); omega_t];
    
    eta(:,i+1) = eta(:,i) + dt * eta_dot(:,i);
    
    cost_value(1,i) = sqrt((x_u0-eta(1,i))^2 + (y_u0-eta(2,i))^2);
end

subplot(2,3,1) % plotting the x position through out the simulation
plot(t, eta(1,1:i),'-r');
set(gca,'fontsize',24)
xlabel('t,[s]');
ylabel('x,[m]');
grid on;

subplot(2,3,2) % plotting the y position through out the simulation
plot(t, eta(2,1:i),'-g');
set(gca,'fontsize',24)
xlabel('t,[s]');
ylabel('y,[m]');
grid on;

subplot(2,3,3)% plotting the shi through out the simulation
plot(t, eta(3,1:i),'-b');
set(gca,'fontsize',24)
xlabel('t,[s]');
ylabel('shi,[rad]');
grid on;

subplot(2,3,4) % plotting the trajectory of target through out the simulation
plot3(eta(1,1:i),eta(2,1:i),z(1:i),'-r');
axis([80 120 140 180 0 200]);
xlabel('X[m]'); ylabel('y[m]'); zlabel('z[m]');
grid on;

subplot(2,3,5) %plotting the cost values
plot(t,cost_value(1,1:i));
xlabel('T[s]');
ylabel('cost value');
grid on;

%% Mobile Robot Animation


l = 0.6; %length of Mobile Robot
w = 0.4; %Width of Mobile Robot
mr_co = [ -l/2, l/2, l/2, -l/2, -l/2;
          -w/2, -w/2, w/2, w/2, -w/2;];
      
subplot(2,3,6)
for i=1:length(t)
    psi = eta(3,i);
    R_psi = [cos(psi), -sin(psi);
             sin(psi), cos(psi);];
    v_pos = R_psi*mr_co;
    
    fill(v_pos(1,:)+eta(1,i), v_pos(2,:)+eta(2,i),'g')
    hold on, grid on
    %axis([50 150 100 200 0 200]), axis square
    plot3(eta(1,1:i),eta(2,1:i),z(1:i),'b-');
    axis([80 120 140 180 0 200]);
    legend('MR','Path')
    set(gca,'fontsize',24)
    xlabel('x,[m]'); ylabel('y[m]');
    zlabel('z,[m]');
    pause(0.1);
    hold off
end

