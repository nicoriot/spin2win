function [ Energy, Mass ] = Energy_Mass_calc(e)
%E_W_CALC System energy and weight calculator for Spin2Win

% Preallocating vectors
Energy = zeros(1,2);
% E(1) = total energy, E(2) = Energy in rpm interval
Mass = zeros(1,e.m+1);
% W(1-m) induvidual shell weight, W(m+1) Total weigth [kg]
I = zeros(1,e.m);
E = zeros(1,m+1);
% E(1->m) = Energy in rpm interval per shell, 
% E(m+1) = Total energy in rpm interval.
% E(m+2) = Total energy in rpm interval per kilogram of rotor.

M = zeros(1,m+1);
% M(1->m) induvidual shell mass, 
% M(m+1) Total mass [kg].

C = zeros(1,m+2);
% C(1->m) = Material Cost of induvidual shell
% C(m+1) = Total material cost
% C(m+2) = Material cost per stored Joule

I = zeros(1,m);

% Useful relations
w = 2*pi*e.n/60;                   % rpm to rad/sec conversion
w_ = 2*pi*e.n_/60;                 % rpm to rad/sec conversion

% Shell mass calculation
for k=1:1:e.m
Mass(k) = e.h*e.p(k)*pi*(e.ro(k)^2-e.ri(k)^2);
end

% Calcualte total mass
Mass(e.m+1) = sum(Mass);

% Shell inertia and energy calculation
for k=1:1:e.m
I(k) = 0.5*Mass(k)*(e.ro(k)^2+e.ri(k)^2);
E(k) = 0.5*I(k)*(w^2-w_^2); %  shell energy in rpm interval
end

% Energy calculation
Energy(m+1) = sum(Energy);
Energy(m+2) = E(m+1)/M(M+1);

% Cost calcualtion
for k=1:1:m
C(k) = M(k)*c(k); %shell cost
end

C(m+1) = sum(C); % total cost
C(m+2) = C(m+1)/E(m+1); % cost per Joule

end

