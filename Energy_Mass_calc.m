function [ Energy, Mass ] = Energy_Mass_calc(e)
%E_W_CALC System energy and weight calculator for Spin2Win

% Preallocating vectors
Energy = zeros(1,2);
% E(1) = total energy, E(2) = Energy in rpm interval
Mass = zeros(1,e.m+1);
% W(1-m) induvidual shell weight, W(m+1) Total weigth [kg]
I = zeros(1,e.m);

% Useful relations
w = 2*pi*e.n/60;                   % rpm to rad/sec conversion
w_ = 2*pi*e.n_/60;                 % rpm to rad/sec conversion

% Shell mass calculation
for k=1:1:e.m
Mass(k) = e.h*e.p(k)*pi*(e.ro(k)^2-e.ri(k)^2);
end

% Calcualte total mass
Mass(e.m+1) = sum(Mass);

% Shell inertia calculation
for k=1:1:e.m
I(k) = 0.5*Mass(k)*(e.ro(k)^2+e.ri(k)^2);
end

% Energy calculation
Energy(1) = 0.5*sum(I)*w^2;
Energy(2) = 0.5*sum(I)*(w^2-w_^2);
end

