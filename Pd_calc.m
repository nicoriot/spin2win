function [z] = Pd_calc(x,e,n)
% Calculate reamaning interferance in geomitry with applied preasssures
% To be used together wiht fsolve and input preasssures data

% Preallocating vectors
dr      = zeros(1,e.m);
Ten_c   = zeros(1,e.m-1);
Ten_r   = zeros(1,e.m-1);
e_z     = zeros(1,e.m);
u       = zeros(1,e.m);
Pi      = zeros(1,e.m-1);
Po      = zeros(1,e.m-1);
C1      = zeros(1,e.m);
C2      = zeros(1,e.m);
Q       = zeros(1,e.m);
Q2      = zeros(1,e.m);
z       = zeros(1,e.m-1);
di      = zeros(1,e.m-1);
do      = zeros(1,e.m-1);
b       = zeros(1,e.m);
a       = zeros(1,e.m);
D       = zeros(1,e.m);
v_rc    = zeros(1,e.m);
v_rz    = zeros(1,e.m);
v_zr    = zeros(1,e.m);
v_zc    = zeros(1,e.m);
v_cz    = zeros(1,e.m);

% Useful relations
w = 2*pi*n/60; % rpm to rad/sec conversion

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
u(k) = sqrt((e.Ec(k)/e.Er(k))*((1-v_rz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k)))); % my
b(k) = (e.v_cr(k)+v_cz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k)); % beta
a(k) = e.Ec(k)*((e_z(k)*(v_zc(k)-v_zr(k)))/(1-v_zc(k)*v_cz(k))); %alpha
end

% Interface calc
for i = 1:1:e.m-1
dr(i) = e.ro(i)-e.ri(i+1);
end
% Make sure interferance is positive
dr = abs(dr);

% Copy input value to preassure vector
Pd = x;

% Format inner and outer preassure vectors using Pd
% Assume same preassure around whole machine
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

D(k) = a(k)/(u(k)^2-1); % D from thesis compliment.

C1(k) = (((e.ri(k)*e.ro(k))^(u(k)))/(e.ri(k)^(2*u(k))-e.ro(k)^(2*u(k))))...
        *(Q(k)*(e.ri(k)^3*e.ro(k)^(u(k))-e.ro(k)^3*e.ri(k)^(u(k)))...
        +(Pi(k)-D(k))*e.ro(k)^(u(k))*e.ri(k)+(D(k)-Po(k))*e.ri(k)^(u(k))*e.ro(k));

C2(k) = (Q(k)*(e.ro(k)^(u(k)+3)-e.ri(k)^(u(k)+3))+(Po(k)-D(k))*e.ro(k)^(u(k)+1)+...
        (D(k)-Pi(k))*e.ri(k)^(u(k)+1))/(e.ri(k)^(2*u(k))-e.ro(k)^(2*u(k)));

end

% Calculate stresses using equations from theory eq 2.20 and 2.22
% And displacements using general Hooke's law
 for k = 1:1:e.m-1
    
    % Calc shell k+1 inner displacement di
    k=k+1;
    
    % Calculate stresses
    Ten_r(k-1) = C1(k)*e.ri(k)^(-1-u(k))+C2(k)*e.ri(k)^(-1+u(k))+Q(k)*e.ri(k)^2-D(k-1);
    
    Ten_c(k-1) = u(k)*(C2(k)*e.ri(k)^(-1+u(k))-C1(k)*e.ri(k)^(-1-u(k)))...
                 +Q2(k)*e.ri(k)^2-D(k-1);
    
    di(k-1) = e.ri(k)*((Ten_c(k-1)/e.Ec(k))*(1-v_zc(k)*v_cz(k))...
              -(v_rc(k)+v_zc(k)*v_rz(k))*Ten_r(k-1)/e.Er(k)-v_zc(k)*e_z(k-1));
          
    k=k-1;
    
    % Calc shell k outer displacement do
    Ten_r(k) = C1(k)*e.ro(k)^(-1-u(k))+C2(k)*e.ro(k)^(-1+u(k))+Q(k)*e.ro(k)^2-D(k);
    
    Ten_c(k) = u(k)*(C2(k)*e.ro(k)^(-1+u(k))-C1(k)*e.ro(k)^(-1-u(k)))...
               +Q2(k)*e.ro(k)^2-D(k);
    
    do(k) = e.ro(k)*((Ten_c(k)/e.Ec(k))*(1-v_zc(k)*v_cz(k))...
            -(v_rc(k)+v_zc(k)*v_rz(k))*Ten_r(k)/e.Er(k)-v_zc(k)*e_z(k)); 
        
 end

% return z as error function
% between inner and outer displacement
% compared to original shell interferance
% scaled up by 10^6 so that fsolve do not end prematurly
 for k=1:1:e.m-1
 z(k) = (di(k)-do(k)-dr(k))*10^6;

 end

end

