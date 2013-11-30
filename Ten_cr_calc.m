function [ Ten_c, Ten_r, intf, e_c, e_r, e_z ] = Ten_cr_calc(x,e,n)
%Calculates hoop and radial stresses of a spinning cylinder

% Ten_c = Hoop stress
% Ten_r = Radial stress
% intf = Radial dicplacement
% e_c = Hoop strain
% e_r = Radial strain
% e_z = Axial strain

% Preallocating vectors
u = zeros(1,e.m);
Pi = zeros(1,e.m-1);
Po = zeros(1,e.m-1);
C1 = zeros(1,e.m);
C2 = zeros(1,e.m);
Q = zeros(1,e.m);
Q2 = zeros(1,e.m);
r = linspace(min(e.ri), max(e.ro), e.y);
Ten_c = zeros(1,length(r));
Ten_r = zeros(1,length(r));
intf = zeros(1,length(r));
e_c = zeros(1,length(r));
e_r = zeros(1,length(r));
e_z = zeros(1,length(r));
b = zeros(1,e.m);
v_rc = zeros(1,e.m);
v_rz = zeros(1,e.m);
v_zr = zeros(1,e.m);
v_zc = zeros(1,e.m);
v_cz = zeros(1,e.m);

% Useful relations
w = 2*pi*n/60;                   % rpm to rad/sec conversion

% Calculate Possions ratios for a assumed transversely isotropic material
for k = 1:1:e.m
v_rc(k) = e.v_cr(k)*e.Er(k)/e.Ec(k);
v_rz(k) = e.Er(k)/(2*e.G_rz(k))-1;
v_zr(k) = v_rz(k);
v_zc(k) = v_rc(k);
v_cz(k) = e.v_cr(k);
end

% Make Young's modulus quotas and b quota
for k = 1:1:e.m
u(k) = sqrt((e.Ec(k)/e.Er(k))*((1-v_rz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k))));
b(k) = (e.v_cr(k)+v_cz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k));
end

% Copy input value to preassure vector
Pd = x;

% Format inner and outer pressure vectors using Pd
% Assume same pressure around whole machine
Pi(1) = 0;
Po(e.m) = 0;
for k=1:1:e.m-1
Pi(k+1) = Pd(k); 
Po(k) = Pd(k);
end

% Make C1 C2 and Q constants using eq 2.28 and 2.29 from theory chapter
for k = 1:1:e.m   
    
Q(k) = (e.p(k)*w^2*(3+b(k))/(u(k)^2-9)); %last part of 2.20
Q2(k) = e.p(k)*w^2*(u(k)^2+3*b(k))/(u(k)^2-9); % last part of 2.22

C1(k) = (((e.ri(k)*e.ro(k))^(u(k)))/(e.ri(k)^(2*u(k))-e.ro(k)^(2*u(k))))...
        *(Q(k)*(e.ri(k)^3*e.ro(k)^(u(k))-e.ro(k)^3*e.ri(k)^(u(k)))...
        +Pi(k)*e.ro(k)^(u(k))*e.ri(k)-Po(k)*e.ri(k)^(u(k))*e.ro(k));

C2(k) = (Q(k)*(e.ro(k)^(u(k)+3)-e.ri(k)^(u(k)+3))+Po(k)*e.ro(k)^(u(k)+1)-...
        Pi(k)*e.ri(k)^(u(k)+1))/(e.ri(k)^(2*u(k))-e.ro(k)^(2*u(k)));

end

% Calculate stresses using equations from theory eq 2.20 and 2.22
% And displacements using general Hooke's law
k = 1;
for i = 1:1:length(r)
    
    k = min(k,e.m-1); %Limit k to prevent index overflow.
    if r(i) > e.ri(k+1)
    k = k+1;
    end
    k = max(k,1); %Limit k to prevent index overflow.
    
    % Calculate stresses
    Ten_r(i) = C1(k)*r(i)^(-1-u(k))+C2(k)*r(i)^(-1+u(k))+Q(k)*r(i)^2;
    
    Ten_c(i) = u(k)*(C2(k)*r(i)^(-1+u(k))-C1(k)*r(i)^(-1-u(k)))+Q2(k)*r(i)^2;
    
    %Calculate displacements
    intf(i) = r(i)*(Ten_c(i)/e.Ec(k)-e.v_cr(k)*Ten_r(i)/e.Er(k));
    % not the correct equations see 2.18 e_z still unknown
end

end

