function [ SE, e_c , e_r, e_z ] = ez_calc(P,e,n)
% strain energy and e_z calcualtor
%
% SE  = Strain Energy
% e_z = axial strain
%
% P = Interface Preassure input data.
% e = data transfer structure.
% n = rotational speed.


disp('*** *** ***    RUNNING EZ_CALC  *** ***  ***');




% Preallocating vectors
Ten_c   = zeros(1,e.m-1);
Ten_r   = zeros(1,e.m-1);
e_z     = zeros(1,e.m);
e_c      = zeros(1,e.m);
e_r      = zeros(1,e.m);
u       = zeros(1,e.m);
Pi      = zeros(1,e.m-1);
Po      = zeros(1,e.m-1);
C1      = zeros(1,e.m);
C2      = zeros(1,e.m);
Q       = zeros(1,e.m);
Q2      = zeros(1,e.m);
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

% Copy input P value to preassure vector
Pd = P;

% Format inner and outer preassure vectors using Pd
% Assume same preassure around whole machine
Pi(1) = 0;
Po(e.m) = 0;
for k=1:1:e.m-1
Pi(k+1) = Pd(k); 
Po(k) = Pd(k);
end

% sweep e_z over range, calcualte energy
e_z =(-0.005:0.00001:0.0025);
SE = zeros(1,length(e_z));
for i = 1:1:length(e_z)
    
% Make Young's modulus quotas and b quota
for k = 1:1:e.m
u(k) = sqrt((e.Ec(k)/e.Er(k))*((1-v_rz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k)))); % my
b(k) = (e.v_cr(k)+v_cz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k)); % beta
a(k) = e.Ec(k)*((e_z(i)*(v_zc(k)-v_zr(k)))/(1-v_zc(k)*v_cz(k))); %alpha
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
 for k = 1:1:e.m
    
    % Calc shell k outer displacement do
    Ten_r(k) = C1(k)*e.ri(k)^(-1-u(k))+C2(k)*e.ri(k)^(-1+u(k))+Q(k)*e.ri(k)^2-D(k);
    
    Ten_c(k) = u(k)*(C2(k)*e.ri(k)^(-1+u(k))-C1(k)*e.ri(k)^(-1-u(k)))...
               +Q2(k)*e.ri(k)^2-D(k);
    
    e_c(i) = (Ten_c(k)/e.Ec(k))*(1-v_zc(k)*v_cz(k))...
            -(Ten_r(k)/e.Er(k))*(v_rc(k)+v_zc(k)*v_rz(k))-v_zc(k)*e_z(i);
          
    e_r(i) = (Ten_r(k)/e.Er(k))*(1-v_zr(k)*v_rz(k))...
            -(Ten_c(k)/e.Ec(k))*(e.v_cr(k)+v_zr(k)*v_cz(k))-v_zr(k)*e_z(i); 
        
 end
 
 % Stress, strain calc done
 % Now calc energy
 
 %SE(i) = (1/(2*(1+v_rz(k))))*(e.Ec(k)*e_c(i)^2+e.Er(k)*e_r(i)^2+...
 %       e.Er(k)*e_z(i)^2+v_rz(k)*(e_c(i)+e_r(i)+e_z(i))/(1-2*v_rz(k))); 
  SE(i) = (1/2)*(e.Ec(k)*e_c(i)^2+e.Er(k)*e_r(i)^2+e.Er(k)*e_z(i)^2)...
          +e.ri(k)^2*w^2*e.p(k)*e_c(i);
           
 
 
end % e_z loop end

end % function end