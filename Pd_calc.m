function [z] = Pd_calc(x,idat,n)
% Calculate reamaning interferance in geomitry with applied preasssures

% Extract useful values from idat
m = idat(8,5);
ri = idat(1,1:1:m);
ro = idat(2,1:1:m);
Ec = idat(3,1:1:m);
Er = idat(4,1:1:m);
v_cr  = idat(5,1:1:m);
p  = idat(6,1:1:m);
G_rz = idat(9,1:1:m);


% Preallocating vectors
dr = zeros(1,m);
Ten_c = zeros(1,m-1);
Ten_r = zeros(1,m-1);
Ten_z = zeros(1,m-1);
u = zeros(1,m);
Pi = zeros(1,m-1);
Po = zeros(1,m-1);
C1 = zeros(1,m);
C2 = zeros(1,m);
Q = zeros(1,m);
Q2 = zeros(1,m);
z = zeros(1,m-1);
di = zeros(1,m-1);
do = zeros(1,m-1);
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

% Interface calc
for i = 1:1:m-1
dr(i) = ro(i)-ri(i+1);
end
% Make sure interferance is positive
dr = abs(dr);

% Copy input value to preassure vector
Pd = x;

% Format inner and outer preassure vectors using Pd
% Assume same preassure around whole machine
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
for k = 1:1:m-1
    
    % Calc shell k+1 inner displacement di
    k=k+1;
    
    Ten_r(k-1) = C1(k)*ri(k)^(-1-u(k))+C2(k)*ri(k)^(-1+u(k))+Q(k)*ri(k)^2;
    
    Ten_c(k-1) = u(k)*(C2(k)*ri(k)^(-1+u(k))-C1(k)*ri(k)^(-1-u(k)))...
                 +Q2(k)*ri(k)^2;
    
    di(k-1) = ri(k)*((Ten_c(k-1)/Ec(k))*(1-v_zc(k)*v_cz(k))...
              -(v_rc(k)+v_zc(k)*v_rz(k))*Ten_r(k-1)/Er(k));
    
    Ten_z(k-1)=Er(k)*(v_zr(k)*Ten_r(k-1)/Er(k)+v_zc(k)*Ten_c(k-1)/Ec(k)); 
          
    k=k-1;
    
    % Calc shell k outer displacement do
    Ten_r(k) = C1(k)*ro(k)^(-1-u(k))+C2(k)*ro(k)^(-1+u(k))+Q(k)*ro(k)^2;
    
    Ten_c(k) = u(k)*(C2(k)*ro(k)^(-1+u(k))-C1(k)*ro(k)^(-1-u(k)))...
               +Q2(k)*ro(k)^2;
    
    do(k) = ro(k)*((Ten_c(k)/Ec(k))*(1-v_zc(k)*v_cz(k))...
            -(v_rc(k)+v_zc(k)*v_rz(k))*Ten_r(k)/Er(k));
    
    Ten_z(k)=Er(k)*(v_zr(k)*Ten_r(k)/Er(k)+v_zc(k)*Ten_c(k)/Ec(k));     
        
end

% return z as error function
% between inner and outer displacement
% compared to original shell interferance
% scaled up so that fsolve do not end prematurly
for k=1:1:m-1
z(k) = (di(k)-do(k)-dr(k))*10^6;
end

end

