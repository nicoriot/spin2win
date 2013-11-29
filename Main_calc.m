function [out] = Main_calc(indata,prefix)
% Pressure and stress calculator for composite and isotropic materials

%%% NOTE!!!!!!!!!!
%This script needs the following m-file functions in order to work:
% Pd_calc.m     - Interface pressure calculator 
% Ten_cr_calc.m - Stress distribution calculator
% E_M_calc      - Energy, Mass calculator
%%%%

% Extract number of active shells from indata
d.m = indata(8,5);

% Place the rest of indata into correct variables
d.ri = indata(1,1:1:d.m)/2000; % [m] Inner radius, scale from mm, dia to m, radius
d.ro = indata(2,1:1:d.m)/2000; % [m] Outer radius, scale from mm, dia to m, radius
d.Ec = indata(3,1:1:d.m)*10^9; % [Pa] Youngs Modulus, Hoop
d.Er = indata(4,1:1:d.m)*10^9; % [Pa ]Youngs Modulus, Radial
d.v  = indata(5,1:1:d.m);      % Possions ratio, cirumference - radial direction
d.p  = indata(6,1:1:d.m);      % [kg/m3] Material density
d.uu = indata(7,1:1:d.m);      % Static friction coifficient
d.G_rz = indata(9,1:1:d.m)*10^9;  % Shear modulus radial-axial
d.C = indata(10,1:1:d.m);  % Shear modulus radial-axial

d.hh = indata(8,2)/1000;     % [m] Cylinder height, scale from mm to m
d.n  = indata(8,3);          % Rotationalspeed rpm
d.n_ = indata(8,4);          % Min speed

% Step size setting, default 0.000001 is enougth for most calcualtions.
step  =  0.000001   ;  % [m] Radial step size used in calculations

% Initiate radial vector to get avaliable flywheel radial positions
r = min(d.ri):step:max(d.ro); % radial vector

% Preallocating vectors
% d.m     = length(d.ri); -maybe not needed due to line 12
dr      = zeros(1,d.m);
Ten_c   = zeros(3,length(r));   
Ten_r   = zeros(3,length(r));
% Where: 1: Pressfit standstill stress, 2: Pure centrifugal stress,
% 3: Pressfit with rotation
intf    = zeros(2,length(r));
% Where: 1: Interface plot standstill, 2: Interface plot rotating
Pd      = zeros(1,d.m-1);            %interface preassure storage vector
K       = zeros(1,d.m-1);
U       = zeros(1,d.m-1);
b       = zeros(1,d.m);
u       = zeros(1,d.m);
j       = zeros(1,d.m);
Ten_o   = zeros(1,d.m);
A       = zeros(3*(d.m-1) , 3*(d.m-1)); 
B       = zeros(3*(d.m-1) , 1);
Pd_still      = zeros(1,d.m-1);
r_mark        = zeros(1,d.m);
r_mark2       = zeros(1,d.m);
warning_flags = zeros(1,5);
Ten_c_max     = zeros(3,d.m);
Ten_c_min     = zeros(3,d.m);
Ten_r_max     = zeros(3,d.m);
Ten_r_min     = zeros(3,d.m);
intf_max      = zeros(2,d.m);
intf_min      = zeros(2,d.m);

e.m     = d.m;
e.ri    = d.ri;
e.ro    = d.ro;
e.Ec    = d.Ec;
e.Er    = d.Er;
e.v_cr  = d.v;
e.p     = d.p;
e.G_rz  = d.G_rz;
e.C     = d.C;
e.h     = d.hh;
e.n     = d.n;
e.n_    = d.n_;
e.y     = length(r);

% Useful constant relations calcualtions
w = 2*pi*d.n/60; % rpm to rad/sec conversion
for k = 1:1:d.m
    % radius quota, used in equation
    b(k) = d.ri(k)/d.ro(k); 

    % µ factor, used in equation
    u(k) = sqrt(d.Ec(k)/d.Er(k));

    % j factor, used in equations
    j(k) = (b(k)^(-u(k)-1)-b(k)^2)/(b(k)^(-u(k)-1)-b(k)^(u(k)-1));

    % Unmodified centrifugal field
    Ten_o(k) = d.p(k)*w^2*d.ro(k)^2*((3+d.v(k))/(9-u(k)^2));
end

% Interface calculation between shells
for i = 1:1:d.m-1
    dr(i) = d.ro(i)-d.ri(i+1);
    if abs(dr(i)) < 0.01*10^-6
        dr(i) = 0.01*10^-6;
        warning_flags(1) = 1;
    end
end

% Make sure interference is positive
dr = abs(dr);

% Calculate stresses due to pure centrifugal forces
% using Ten_cr_calc subroutine with zero shell interference
Pdz = zeros(1,d.m-1);

[Ten_c(2,1:1:length(r)), Ten_r(2,1:1:length(r)), ...
intf(2,1:1:length(r))] = Ten_cr_calc(Pdz,e,d.n);

%skip all pressfit calcualtions if there is less than 2 shells
if d.m ~= 1

% Calculate approximate preassures using 
% a simple linear matrix. Based on standard pressfit formulas.

% Section preassure calc factors used i matrix
for k = 1:1:d.m-1
    K(k) = d.ri(k+1)*((d.ro(k+1)^2+d.ri(k+1)^2)/(d.Ec(k+1)*(d.ro(k+1)^2 ...
          -d.ri(k+1)^2))+d.v(k+1)/d.Er(k+1));
    U(k) = d.ro(k)*((d.ro(k)^2+d.ri(k)^2)/(d.Ec(k)*(d.ro(k)^2-d.ri(k)^2))-d.v(k)/d.Er(k));    
end

% Example Matrix structure for a 4 shelled geometry.

% B = [ 0      ;
%       0      ; 
%       0      ; 
%       0      ; 
%       0      ; 
%       0      ; 
%       dr(1)  ; 
%       dr(2)  ;
%       dr(3)  ];

% %     P    P    P     o     o    o    i    i    i
% A = [ K(1) 0    0    -1     0    0    0    0    0  ;
%       0    K(2) 0     0    -1    0    0    0    0  ;
%       0    0    K(3)  0     0   -1    0    0    0  ;
%      -U(1) 0    0     0     0    0   -1    0    0  ;
%       0   -U(2) 0     0     0    0    0   -1    0  ;
%       0    0   -U(3)  0     0    0    0    0   -1  ;
%       0    0    0     1     0    0   -1    1    0  ;
%       0    0    0    -1     1    0    0   -1    1  ;
%       0    0    0     0    -1    1    0    0   -1  ];  


% Matrix maker

% Pd do fill
for k=1:1:d.m-1
    A(k , k)=K(k);
    A(k , (d.m-1+k)) = -1;  
end

% Pd di fill
for k=1:1:d.m-1
    A(d.m-1+k , k) = -U(k);
    A(d.m-1+k , 2*(d.m-1)+k) = -1;
end

% do di fill 1
for k=1:1:d.m-1
    A(2*(d.m-1)+k , d.m-1+k) = 1;
    A(2*(d.m-1)+k , 2*(d.m-1)+k) = -1;  
end

% do di fill 2
for k=1:1:d.m-2
    A(2*(d.m-1)+1+k , d.m-1+k) = -1;
    A(2*(d.m-1)+k , 2*(d.m-1)+1+k) = 1;
end

% B fill
for k=1:1:d.m-1
    B(2*(d.m-1)+k) = dr(k);  
end

% Matrix solver
D = A\B;

% Extract preassures from vector X
% Pd will be approximate now if a isotropic material is used.
% Values will be used as input guess into fsolver for finetuning 
%or major change if an orthotropic material is used  
Pd_t = 0;
for k = 1:1:d.m-1
    Pd_t(k) = D(k);
end
  
% Copy values to pressure vector
Pd = Pd_t;

% Special Fsolve settings for increased chance of convergence 
options = optimset('TolFun',1e-9,'TolFun',1e-9);

% Finetune pressure guess from above matrix calculations
% by using fsolve on Pd_calc subroutine
% n = 0 for standstill curve
[Pd, ~,warning_flags(4)] = fsolve(@(x) Pd_calc(x,e,0),Pd,options);
Pd_still = Pd;

% Calcualte the stresses in standstill case
% with the Ten_cr_calc subroutine
[Ten_c(1,1:1:length(r)), Ten_r(1,1:1:length(r)), ...
    intf(1,1:1:length(r))] = Ten_cr_calc(Pd,e,0);

% Reset Pd to previous approximate values for next run
Pd = Pd_t;

% Finetune preassure guess, include centrifugal pressure, n =/= 0
[Pd, ~,warning_flags(5)] = fsolve(@(x) Pd_calc(x,e,d.n),Pd,options);

% Calcualte Pressfit stresses due to rotating pressfit
% with the Ten_cr_calc subroutine
[Ten_c(3,1:1:length(r)), Ten_r(3,1:1:length(r)), ...
     intf(2,1:1:length(r))] = Ten_cr_calc(Pd,e,d.n);
end % end for if statement about number of shells

% Sacling from Pa to MPa for easier enterpritation
Ten_c = Ten_c/(10^6);
Ten_r = Ten_r/(10^6);

% Find shell surface radius r length value
% so that shell interfaces can be marked

for i=1:1:d.m-1
r_mark(1,i) = (d.ri(i+1)/step-d.ri(1)/step)+1;

end
r_mark(1,d.m)=(d.ro(d.m)/step-d.ri(1)/step)+1;
r_mark = floor(r_mark);

for i=1:1:d.m-1
    r_mark2(1,i+1)=r_mark(1,i)+2;
end
r_mark2(1,1)=1;

% Make min max values vectors for each r_mark
for i=1:1:3
    for k=1:1:d.m
      Ten_c_max(i,k) = max(Ten_c(i,(r_mark2(k)+2):r_mark(k)));
      Ten_c_min(i,k) = min(Ten_c(i,(r_mark2(k)+2):r_mark(k)));
    end
end

for i=1:1:3
    for k=1:1:d.m
      Ten_r_max(i,k) = max(Ten_r(i,(r_mark2(k)+2):r_mark(k)));
      Ten_r_min(i,k) = min(Ten_r(i,(r_mark2(k)+2):r_mark(k)));
    end
end

for i=1:1:2
    for k=1:1:d.m
      intf_max(i,k) = max(intf(i,(r_mark2(k)+2):r_mark(k)));
      intf_min(i,k) = min(intf(i,(r_mark2(k)+2):r_mark(k)));
    end
end

% Calculate mass and energy of system using E_M_calc subroutine
[Energy, Mass, Cost] = Energy_Mass_calc(e);


% Print data to screen and save into logfile named after the current date
% Mark each run with a timestamp

% Start diary mode
diary(['Spin2Win_output_',date,'.txt'])
diary on

% Get date and time from system
stamp = datestr(now,0);
disp(['Date and time of run: ', stamp, ', ', prefix])

% Print input data
format long
disp('Input data:')
disp(['Ri: ', num2str(d.ri,10)])
disp(['Ro: ', num2str(d.ro,10)])
disp(['Ec: ', num2str(d.Ec,10)])
disp(['Er: ', num2str(d.Er,10)])
disp(['v: ', num2str(d.v,10)])
disp(['p: ', num2str(d.p,10)])
disp(['µ: ', num2str(d.uu,10)])
disp(['G: ', num2str(d.G_rz,10)])
disp(['C: ', num2str(d.C,10)])
disp(['lr: ', num2str(length(r),10)])
disp(['h: ', num2str(d.hh,10)])
disp(['n: ', num2str(d.n,10)])
disp(['n_: ', num2str(d.n_,10)])
disp(['sh: ', num2str(d.m,10)])
fprintf('\n')
format short

% Print output data
disp('Results:')
disp('--------------------------------------------------')
disp('Centrifugal only:')
disp(['Maximal tension, hoop: ', num2str(max(Ten_c(2,:))), ' [MPa]'])
disp(['Minimal tension, hoop: ', num2str(min(Ten_c(2,:))), ' [MPa]'])
disp(['Maximal tension, radial: ', num2str(max(Ten_r(2,:))), ' [MPa]'])
disp(['Minimal tension, radial: ', num2str(min(Ten_r(2,:))), ' [MPa]'])
disp('--------------------------------------------------')
disp('Press-fit, standstill:')
disp(['Maximal tension, hoop: ', num2str(max(Ten_c(1,:))), ' [MPa]'])
disp(['Minimal tension, hoop: ', num2str(min(Ten_c(1,:))), ' [MPa]'])
disp(['Maximal tension, radial: ', num2str(max(Ten_r(1,:))), ' [MPa]'])
disp(['Minimal tension, radial: ', num2str(min(Ten_r(1,:))), ' [MPa]'])
disp('Interference pressures:')
%Repeat for number of interferences
for k = 1:1:d.m-1
   disp(['Shell ', num2str(k), '-' , num2str(k+1), ': ', ...
       num2str(Pd_still(k)/10^6), ' [MPa]']) 
end
disp('--------------------------------------------------')
disp('Press-fit, rotating:')
disp(['Maximal tension, hoop: ', num2str(max(Ten_c(3,:))), ' [MPa]'])
disp(['Minimal tension, hoop: ', num2str(min(Ten_c(3,:))), ' [MPa]'])
disp(['Maximal tension, radial: ', num2str(max(Ten_r(3,:))), ' [MPa]'])
disp(['Minimal tension, radial: ', num2str(min(Ten_r(3,:))), ' [MPa]'])
disp('Interference pressures:')
%Repeat for number of interferences
for k = 1:1:d.m-1
   disp(['Shell ', num2str(k), '-' , num2str(k+1), ': ', ... 
       num2str(Pd(k)/10^6), ' [MPa]']) 
end
disp('--------------------------------------------------')
disp('Maximal assembly force needed for shell press-fit mounting:')
%Repeat for number of interferences
for k = 1:1:d.m-1
    F = pi*2*d.ro(k)*d.hh*d.uu(k)*Pd_still(k)/1000;
    disp(['Shell ', num2str(k), '-' , num2str(k+1), ': ', num2str(F), ' [kN]']) 
end
disp('--------------------------------------------------')
disp('Maximal transferable torque between shells, at standstill:')
%Repeat for number of interferences
for k = 1:1:d.m-1
    F = pi*2*d.ro(k)*d.hh*d.uu(k)*Pd_still(k)*d.ro(k)/1000;
    disp(['Shell ', num2str(k), '-' , num2str(k+1), ': ', num2str(F), ' [kNm]']) 
end
disp('--------------------------------------------------')
disp('Maximal transferable torque between shells, at max rpm:')
%Repeat for number of interferences
for k = 1:1:d.m-1
    F = pi*2*d.ro(k)*d.hh*d.uu(k)*Pd(k)*d.ro(k)/1000;
    disp(['Shell ', num2str(k), '-' , num2str(k+1), ': ', num2str(F), ' [kNm]']) 
end
disp('--------------------------------------------------')

disp('Shell and system MASS:')
%Repeat for number of shells
for k = 1:1:d.m
    disp(['Shell ', num2str(k), ': ', num2str(Mass(k)), ' [kg]']) 
end
disp(['Total: ', num2str(Mass(d.m+1)), ' [kg]']) 
disp('--------------------------------------------------')

disp('Shell and system ENERGY:')
%Repeat for number of shells
for k = 1:1:d.m
    disp(['Shell ', num2str(k), ': ', num2str(Energy(k)/3600), ' [Wh]']) 
end
disp(['Total: ', num2str(Energy(d.m+1)/3600), ' [Wh]']) 
disp(['Total per kg: ', num2str(Energy(d.m+2)/3600), ' [Wh/kg]'])
disp('--------------------------------------------------')

disp('Shell and system COST:')
%Repeat for number of shells
for k = 1:1:d.m
    disp(['Shell ', num2str(k), ': ', num2str(Cost(k)), ' [$]']) 
end
disp(['Total: ', num2str(Cost(d.m+1)), ' [$]']) 
disp(['Total per Wh: ', num2str(Cost(d.m+2)*3600), ' [$/Wh]'])
disp('--------------------------------------------------')


% Print WARNINGS if they have been flagged
if min(Pd) < 0
    warning_flags(2) = 1;
end
if min(Pd_still) < 0
    warning_flags(3) = 1;
end

if warning_flags(1) == 1
disp('*************WARNING******WARNING**************')
disp('Warning: Interference non existent or to small. < 0.1*10^-6')
disp('Emergency minimum value of: 0.1*10^-6 [m] activated.')
disp('*************WARNING******WARNING**************')
end

if warning_flags(2) == 1
disp('*************WARNING******WARNING**************')
disp('Warning: Rotating interference pressure have gone negative,')
disp('calculations higly unstable and mostly broken.')
disp('Reduce speed or increase interference')
disp('*************WARNING******WARNING**************')
end

if warning_flags(3) == 1
disp('*************WARNING******WARNING**************')
disp('Warning: Standstill interference pressure have gone negative,')
disp('calculations higly unstable and mostly broken.')
disp('Check input data for sign errors and interference')
disp('*************WARNING******WARNING**************')
end

if warning_flags(4) ~= 1
disp('*************WARNING******WARNING**************')
disp('Warning: Standstill pressure calculation did not converge,')
disp('Check input data for sign errors and interference')
disp('*************WARNING******WARNING**************')
end

if warning_flags(5) ~= 1
disp('*************WARNING******WARNING**************')
disp('Warning: Rotating pressure calculation did not converge,')
disp('Check input data for sign errors and interference')
disp('*************WARNING******WARNING**************')
end

% Mark end of run and fill wiht some blank lines.
disp('End of output')
fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('\n')
fprintf('\n')

% Stop diary mode
diary off


% Prepere and Copy important data and plotter data into struct for return
% Plotter data
out.Ten_c = Ten_c     ;
out.Ten_r = Ten_r     ;
out.intf = intf       ;
out.r_mark = r_mark   ;
out.r_mark2 = r_mark2 ;
out.r = r             ;
out.ri = d.ri           ;
out.stamp = stamp     ;

% General data
out.P_standstill = Pd_still  ;
out.P_rotating   = Pd        ;
out.Energy       = Energy    ;
out.Mass         = Mass      ;
out.Ten_c_max    = Ten_c_max ;
out.Ten_c_min    = Ten_c_min ;
out.Ten_r_max    = Ten_r_max ;
out.Ten_r_min    = Ten_r_min ;

end
