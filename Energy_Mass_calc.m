function [ E, M, C] = Energy_Mass_calc(e)
%E_W_CALC System energy and weight calculator for Spin2Win

% Preallocating vectors
I = zeros(1,e.m);
E = zeros(1,e.m+2);
% E(1->m) = Energy in rpm interval per shell, 
% E(m+1) = Total energy in rpm interval.
% E(m+2) = Total energy in rpm interval per kilogram of rotor.

M = zeros(1,e.m+1);
% M(1->m) induvidual shell mass, 
% M(m+1) Total mass [kg].

C = zeros(1,e.m+2);
% C(1->m) = Material Cost of induvidual shell
% C(m+1) = Total material cost
% C(m+2) = Material cost per stored Joule

% Useful relations
w = 2*pi*e.n/60;                   % rpm to rad/sec conversion
w_ = 2*pi*e.n_/60;                 % rpm to rad/sec conversion

% Shell mass calculation
for k=1:1:e.m
M(k) = e.h*e.p(k)*pi*(e.ro(k)^2-e.ri(k)^2);
end

% Calcualte total mass
M(e.m+1) = sum(M);

% Shell inertia and energy calculation
for k=1:1:e.m
I(k) = 0.5*M(k)*(e.ro(k)^2+e.ri(k)^2);
E(k) = 0.5*I(k)*(w^2-w_^2); %  shell energy in rpm interval
end

% Energy calculation
E(e.m+1) = sum(E);
E(e.m+2) = E(e.m+1)/M(e.m+1);

% Cost calcualtion
for k=1:1:e.m
C(k) = M(k)*e.C(k); %shell cost
end

C(e.m+1) = sum(C); % total cost
C(e.m+2) = C(e.m+1)/E(e.m+1); % cost per Joule

end

