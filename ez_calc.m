function [ SE, e_zmin ] = ez_calc(P,e,n,g)
% strain energy minimizer and e_z calcualtor
%
% SE  = Strain Energy
% e_zmin = axial strain at min energy
%
% P = Interface Preassure input data.
% e = data transfer structure.
% n = rotational speed.
% g = current k (shell) of parent loop

disp('*** *** ***  RUNNING EZ_CALC  *** *** ***');

% Preallocating vectors
Ten_c   = zeros(1,e.m-1);
Ten_r   = zeros(1,e.m-1);
Ten_z   = zeros(1,e.m-1);
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

k = g; %copy g as current k (shell number)

% Calculate Possions ratios for a assumed transversely isotropic material
v_rc(k) = e.v_cr(k)*e.Er(k)/e.Ec(k);
v_rz(k) = e.Er(k)/(2*e.G_rz(k))-1;
v_zr(k) = v_rz(k);
v_zc(k) = v_rc(k);
v_cz(k) = e.v_cr(k);

% sweep e_z over range, calcualte energy
e_z      =(-0.005:0.00001:-0.001122);
SE       = zeros(1,length(e_z));
e_c      = zeros(1,length(e_z));
e_r      = zeros(1,length(e_z));


for i = 1:1:length(e_z) % e_z sweep loop
    
% note: k still used as index despite that e_z calc only use one k at a
% time. This is a prepiration for if many k wants to be run in the samme
% loop.
    
% Make Young's modulus quotas, a and b quota
u(k) = sqrt((e.Ec(k)/e.Er(k))*((1-v_rz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k)))); % my
b(k) = (e.v_cr(k)+v_cz(k)*v_zr(k))/(1-v_zc(k)*v_cz(k)); % beta
a(k) = e.Ec(k)*((e_z(i)*(v_zc(k)-v_zr(k)))/(1-v_zc(k)*v_cz(k))); %alpha


% Make C1 C2 and Q constants using eq 2.28 and 2.29 from theory chapter 
% Complimented C1 and C2 from thesis compliment    
Q(k) = (e.p(k)*w^2*(3+b(k))/(u(k)^2-9)); %last part of 2.20
Q2(k) = e.p(k)*w^2*(u(k)^2+3*b(k))/(u(k)^2-9); % last part of 2.22
D(k) = a(k)/(u(k)^2-1); % D from thesis compliment.

C1(k) = (((e.ri(k)*e.ro(k))^(u(k)))/(e.ri(k)^(2*u(k))-e.ro(k)^(2*u(k))))...
        *(Q(k)*(e.ri(k)^3*e.ro(k)^(u(k))-e.ro(k)^3*e.ri(k)^(u(k)))...
        +(Pi(k)-D(k))*e.ro(k)^(u(k))*e.ri(k)+(D(k)-Po(k))*e.ri(k)^(u(k))*e.ro(k));

C2(k) = (Q(k)*(e.ro(k)^(u(k)+3)-e.ri(k)^(u(k)+3))+(Po(k)-D(k))*e.ro(k)^(u(k)+1)+...
        (D(k)-Pi(k))*e.ri(k)^(u(k)+1))/(e.ri(k)^(2*u(k))-e.ro(k)^(2*u(k)));

    
% Calc shell stress and strains
Ten_r(k) = C1(k)*e.ri(k)^(-1-u(k))+C2(k)*e.ri(k)^(-1+u(k))+Q(k)*e.ri(k)^2-D(k);
    
Ten_c(k) = u(k)*(C2(k)*e.ri(k)^(-1+u(k))-C1(k)*e.ri(k)^(-1-u(k)))...
           +Q2(k)*e.ri(k)^2-D(k);

Ten_z(k) = e.Er(k)*(e_z(i)+v_rz(k)*Ten_r(k)/e.Er(k)+v_cz(k)*Ten_c(k)/e.Ec(k));
    
e_c(i) = (Ten_c(k)/e.Ec(k))*(1-v_zc(k)*v_cz(k))...
         -(Ten_r(k)/e.Er(k))*(v_rc(k)+v_zc(k)*v_rz(k))-v_zc(k)*e_z(i);
          
e_r(i) = (Ten_r(k)/e.Er(k))*(1-v_zr(k)*v_rz(k))...
         -(Ten_c(k)/e.Ec(k))*(e.v_cr(k)+v_zr(k)*v_cz(k))-v_zr(k)*e_z(i);

u_r(k) = e_c(i)*e.ri(k);
u_c(k) = 0;
u_z(k) = 0;
        
 % Stress, strain calc done
 
 % Calc energy for current strain. Correct????
 % SE(i) = (1/2)*(e.Ec(k)*e_c(i)^2+e.Er(k)*e_r(i)^2+e.Er(k)*e_z(i)^2)...
 %         +e.ri(k)^2*w^2*e.p(k)*e_c(i);
 
  SE(i) = (1/2)*(Ten_c(k)*e_c(i)+Ten_r(k)*e_r(i)+Ten_z(k)*e_z(i))...
          -e.ri(k)^2*w^2*e.p(k)*e_c(i);
     
 % SE(i) = 5.59e5...
 %       -e.ri(k)^2*w^2*e.p(k)*e_c(i);
     
   
end % e_z loop end
  
  %%% HAXXPLOT rememeber to remove
  figure(5)
  plot(e_z,SE,'-b')
  hold on
  
  
 TEN_CC = Ten_c(k)
 TEN_RR = Ten_r(k)
 TEN_ZZ = Ten_z(k)
 EPS_CC = e_c(i)
 EPS_RR = e_r(i)
 EPS_ZZ = e_z(i)
 RADIUS = e.ri(k)
 RAD_SPEED = w
 DENSITY = e.p(k)
 STRAIN_E_DENS = 0.5*(Ten_c(k)*e_c(i)+Ten_r(k)*e_r(i)+Ten_z(k)*e_z(i))
 EXTERNAL_E_DENS = e.ri(k)^2*w^2*e.p(k)*e_c(i)
 
 %%% HAXXPLOT

  [~,N] = min(SE); % find min energy in interval
  e_zmin = -e_z(N); % what is min energy's e_z, 
  plot(e_z(N),SE(N),'*r')
  
  %WHY -e_z????? It works for iso materials but breaks for ortho.
  %Check energy formula.

end % function end