% Target state pridiction model
clear all; clc; close all;

%Simulation Parameters
dt = 0.1; %tim-size
ts = 10; %Simulation time
t = 0:dt:ts; %Time span

%initial location of the target according to paper
x0 = 100;
y0 = 150;
shi0 = 0;

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
    
end

figure % plotting the x position through out the simulation
plot(t, eta(1,1:i),'-r');
set(gca,'fontsize',24)
xlabel('t,[s]');
ylabel('x,[m]');
grid on;

figure % plotting the y position through out the simulation
plot(t, eta(2,1:i),'-g');
set(gca,'fontsize',24)
xlabel('t,[s]');
ylabel('y,[m]');
grid on;

figure % plotting the shi through out the simulation
plot(t, eta(3,1:i),'-b');
set(gca,'fontsize',24)
xlabel('t,[s]');
ylabel('shi,[rad]');
grid on;

figure % plotting the trajectory of target through out the simulation
plot3(eta(1,1:i),eta(2,1:i),z(1:i),'-r');
grid on;

