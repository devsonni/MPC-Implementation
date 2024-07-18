function [x,y] = ellipse(a_p, b_p, X_E_1, Y_E_1)
% hold on
x = [];
y = [];
th = 0:pi/50:2*pi;
x = a_p*sin(th)+X_E_1;
y = b_p*cos(th)+Y_E_1;
% hold off
end
