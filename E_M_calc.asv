function [ E, M ] = E_M_calc( idat )
%E_W_CALC System energy and weight calculator for Spin2Win

% Extract values from idat
m  = idat(8,5);
ri = (idat(1,1:1:m));
ro = (idat(2,1:1:m));
p  = idat(6,1:1:m);
h  = (idat(8,2));
n  = idat(8,3);
n_ = idat(8,4);

% Preallocating vectors
E = zeros(1,2);
% E(1) = total energy, E(2) = Energy in rpm interval
M = zeros(1,m+1);
% W(1-m) induvidual shell weight, W(m+1) Total weigth [kg]
I = zeros(1,m);

% Useful relations
w = 2*pi*n/60;                   % rpm to rad/sec conversion
w_ = 2*pi*n_/60;                 % rpm to rad/sec conversion

% Shell mass calculation
for k=1:1:m
M(k) = h*p(k)*pi*(ro(k)^2-ri(k)^2);
end

% Calcualte total mass
M(m+1) = sum(M);

% Shell inertia calculation
for k=1:1:m
I(k) = 0.5*M(k)*(ro(k)^2+ri(k)^2);
end

% Energy calculation
E(1) = 0.5*sum(I)*w^2;
E(2) = 0.5*sum(I)*(w^2-w_^2);
end

