function [ E, M ] = E_M_calc(e)
%E_W_CALC System energy and weight calculator for Spin2Win

% Extract values from idat
% e.m  = idat(8,5);
% e.ri = (idat(1,1:1:e.m));
% e.ro = (idat(2,1:1:e.m));
% e.p  = idat(6,1:1:e.m);
% e.h  = (idat(8,2));
% e.n  = idat(8,3);
% e.n_ = idat(8,4);

% Preallocating vectors
E = zeros(1,2);
% E(1) = total energy, E(2) = Energy in rpm interval
M = zeros(1,e.m+1);
% W(1-m) induvidual shell weight, W(m+1) Total weigth [kg]
I = zeros(1,e.m);

% Useful relations
w = 2*pi*e.n/60;                   % rpm to rad/sec conversion
w_ = 2*pi*e.n_/60;                 % rpm to rad/sec conversion

% Shell mass calculation
for k=1:1:e.m
M(k) = e.h*e.p(k)*pi*(e.ro(k)^2-e.ri(k)^2);
end

% Calcualte total mass
M(e.m+1) = sum(M);

% Shell inertia calculation
for k=1:1:e.m
I(k) = 0.5*M(k)*(e.ro(k)^2+e.ri(k)^2);
end

% Energy calculation
E(1) = 0.5*sum(I)*w^2;
E(2) = 0.5*sum(I)*(w^2-w_^2);
end

