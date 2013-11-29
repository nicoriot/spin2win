function [ Ten_c, Ten_r, intf ] = Ten_cr_calc(x,idat,n)
%Calculates hoop and radial stresses of a spinning cylinder

% Extract values from idat
m = idat(8,5);
ri = idat(1,1:1:m);
ro = idat(2,1:1:m);
Ec = idat(3,1:1:m);
Er = idat(4,1:1:m);
v_cr  = idat(5,1:1:m);
p  = idat(6,1:1:m);
y  = idat(8,1);
G_rz = idat(9,1:1:m);

% Preallocating vectors
u = zeros(1,m);
Pi = zeros(1,m-1);
Po = zeros(1,m-1);
C1 = zeros(1,m);
C2 = zeros(1,m);
Q = zeros(1,m);
Q2 = zeros(1,m);
r = linspace(min(ri), max(ro), y);
Ten_c = zeros(1,length(r));
Ten_r = zeros(1,length(r));
intf = zeros(1,length(r));
b = zeros(1,m);
v_rc = zeros(1,m);
v_rz = zeros(1,m);
v_zr = zeros(1,m);
v_zc = zeros(1,m);
v_cz = zeros(1,m);

% Useful relations
w = 2*pi*n/60;                   % rpm to rad/sec conversion

% Calculate Possions ratios for a assumed transversely isotropic material
for k = 1:1:m
v_rc(k) = v_cr(k)*Er(k)/Ec(k);
v_rz(k) = Er(k)/(2*G_rz(k))-1;
v_zr(k) = v_rz(k);
v_zc(k) = v_rc(k);
v_cz(k) = v_cr(k);
end

% Make Young's modulus quotas and b quota
for k = 1:1:m
u(k) = sqrt((Ec(k)/Er(k))*((1-v_rz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k))));
b(k) = (v_cr(k)+v_cz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k));
end

% Copy input value to preassure vector
Pd = x;

% Format inner and outer pressure vectors using Pd
% Assume same pressure around whole machine
Pi(1) = 0;
Po(m) = 0;
for k=1:1:m-1
Pi(k+1) = Pd(k); 
Po(k) = Pd(k);
end

% Make C1 C2 and Q constants using eq 2.28 and 2.29 from theory chapter
for k = 1:1:m   
    
Q(k) = (p(k)*w^2*(3+b(k))/(u(k)^2-9)); %last part of 2.20
Q2(k) = p(k)*w^2*(u(k)^2+3*b(k))/(u(k)^2-9); % last part of 2.22

C1(k) = (((ri(k)*ro(k))^(u(k)))/(ri(k)^(2*u(k))-ro(k)^(2*u(k))))...
        *(Q(k)*(ri(k)^3*ro(k)^(u(k))-ro(k)^3*ri(k)^(u(k)))...
        +Pi(k)*ro(k)^(u(k))*ri(k)-Po(k)*ri(k)^(u(k))*ro(k));

C2(k) = (Q(k)*(ro(k)^(u(k)+3)-ri(k)^(u(k)+3))+Po(k)*ro(k)^(u(k)+1)-...
        Pi(k)*ri(k)^(u(k)+1))/(ri(k)^(2*u(k))-ro(k)^(2*u(k)));

end

% Calculate stresses using equations from theory eq 2.20 and 2.22
% And displacements using general Hooke's law
k = 1;
for i = 1:1:length(r)
    
    k = min(k,m-1); %Limit k to prevent index overflow.
    if r(i) > ri(k+1)
    k = k+1;
    end
    k = max(k,1); %Limit k to prevent index overflow.
    
    % Calculate stresses
    Ten_r(i) = C1(k)*r(i)^(-1-u(k))+C2(k)*r(i)^(-1+u(k))+Q(k)*r(i)^2;
    
    Ten_c(i) = u(k)*(C2(k)*r(i)^(-1+u(k))-C1(k)*r(i)^(-1-u(k)))+Q2(k)*r(i)^2;
    
    %Calculate displacements
    intf(i) = r(i)*(Ten_c(i)/Ec(k)-v_cr(k)*Ten_r(i)/Er(k));
end

end

