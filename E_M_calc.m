function [ E, M, C ] = E_M_Calc( idat )
%E_W_CALC System energy, weight and cost calculator for Spin2Win

%Calcualtions made in Joule, convert in Main_calc to Wh if wanted.

% Extract values from idat
m  = idat(8,5);
ri = idat(1,1:1:m);
ro = idat(2,1:1:m);
p  = idat(6,1:1:m);
c  = idat(10,1:1:m);
h  = idat(8,2);
n  = idat(8,3);
n_ = idat(8,4);

% Preallocating vectors
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
w = 2*pi*n/60;                   % rpm to rad/sec conversion
w_ = 2*pi*n_/60;                 % rpm to rad/sec conversion

% Shell mass calculation
for k=1:1:m
M(k) = h*p(k)*pi*(ro(k)^2-ri(k)^2);
end

% Calcualte total mass
M(m+1) = sum(M);

% Shell inertia and energy calculation
for k=1:1:m
I(k) = 0.5*M(k)*(ro(k)^2+ri(k)^2); % inertia
E(k) = 0.5*I(k)*(w^2-w_^2); %  shell energy in rpm interval
end

% Energy calculation
E(m+1) = sum(E); % energy in rpm interval
E(m+2) = E(m+1)/M(m+1);

% Cost calcualtion
for k=1:1:m
C(k) = M(k)*c(k); %shell cost
end

C(m+1) = sum(C); % total cost
C(m+2) = C(m+1)/E(m+1); % cost per Joule

end

