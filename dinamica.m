function dxdt = dinamica(t,x)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%G = -9.8;
%B = 0.1 ;
%J = 0.01;
%L = 0.005;
%R = 0.01;
%K = 0.01;
u = 1;
%r = 2;

G = 9.8;
B = 0.01;
J = 3;
L = 0.001;
R = 0.001;
K = 0.0879;
r = 2;

dxdt = [x(2) ; sin(x(1)) * G / (J * r)   - x(2)* B / (J * r) + K * x(3); - x(2)* K/ L - x(3) * R / L + u / L];
end

