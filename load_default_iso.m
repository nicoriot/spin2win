function [ idat,Data ] = load_default_iso()
% Syntax:
% [Output data] = load_default_iso()

% Generates a set of data that can be used as input data in Spin2Win

% Input paramaters
%                  Aluminum SETTINGS
%       C1          C2          C3          C4         C5 
%       C6          C7          C8          C9 
% [m] Inner diameter
r_scale = 2000;
Data.Shells(1).d.i = 0.100*r_scale;
Data.Shells(2).d.i = 0.125*r_scale;
Data.Shells(3).d.i = 0.150*r_scale;
Data.Shells(4).d.i = 0.175*r_scale;
Data.Shells(5).d.i = 0.200*r_scale;
Data.Shells(6).d.i = 0.225*r_scale;
Data.Shells(7).d.i = 0.250*r_scale;
Data.Shells(8).d.i = 0.275*r_scale;
Data.Shells(9).d.i = 0.300*r_scale;

ri = [ 0.100    ,  0.125    ,  0.150    ,  0.175    ,  0.200     , ...  
       0.225    ,  0.250    ,  0.275    ,  0.300    ]; 

% [m] Outer diameter
Data.Shells(1).d.o = 0.12525*r_scale;
Data.Shells(2).d.o = 0.15025*r_scale;
Data.Shells(3).d.o = 0.17525*r_scale;
Data.Shells(4).d.o = 0.20025*r_scale;
Data.Shells(5).d.o = 0.22525*r_scale;
Data.Shells(6).d.o = 0.25025*r_scale;
Data.Shells(7).d.o = 0.27525*r_scale;
Data.Shells(8).d.o = 0.30025*r_scale;
Data.Shells(9).d.o = 0.32525*r_scale;

ro = [ 0.12525  ,  0.15025  ,  0.17525  ,  0.20025  ,  0.22525   , ... 
       0.25025  ,  0.27525  ,  0.30025  ,  0.32500  ]; 
   
% Youngs Modulus, Hoop [MPa]
Data.Shells(1).E.c = 72;
Data.Shells(2).E.c = 72;
Data.Shells(3).E.c = 72;
Data.Shells(4).E.c = 72;
Data.Shells(5).E.c = 72;
Data.Shells(6).E.c = 72;
Data.Shells(7).E.c = 72;
Data.Shells(8).E.c = 72;
Data.Shells(9).E.c = 72;

Ec = [ 72*10^9  ,  72*10^9  ,  72*10^9  ,  72*10^9  ,  72*10^9   , ... 
       72*10^9  ,  72*10^9  ,  72*10^9  ,  72*10^9  ]; 

% Youngs Modulus, Radial [MPa] 
Data.Shells(1).E.r = 72;
Data.Shells(2).E.r = 72;
Data.Shells(3).E.r = 72;
Data.Shells(4).E.r = 72;
Data.Shells(5).E.r = 72;
Data.Shells(6).E.r = 72;
Data.Shells(7).E.r = 72;
Data.Shells(8).E.r = 72;
Data.Shells(9).E.r = 72;

Er = [ 72*10^9  ,  72*10^9  ,  72*10^9  ,  72*10^9  ,  72*10^9   , ... 
       72*10^9  ,  72*10^9  ,  72*10^9  ,  72*10^9  ]; 
% Possions ratio, cirumference - radial direction   
Data.Shells(1).E.v = 0.33;
Data.Shells(2).E.v = 0.33;
Data.Shells(3).E.v = 0.33;
Data.Shells(4).E.v = 0.33;
Data.Shells(5).E.v = 0.33;
Data.Shells(6).E.v = 0.33;
Data.Shells(7).E.v = 0.33;
Data.Shells(8).E.v = 0.33;
Data.Shells(9).E.v = 0.33;

v  = [ 0.33     ,  0.33     ,  0.33     ,  0.33     ,  0.33      , ... 
       0.33     ,  0.33     ,  0.33     ,  0.33     ]; 

% [kg/m3] Material density   
Data.Shells(1).rho = 2800;
Data.Shells(2).rho = 2800;
Data.Shells(3).rho = 2800;
Data.Shells(4).rho = 2800;
Data.Shells(5).rho = 2800;
Data.Shells(6).rho = 2800;
Data.Shells(7).rho = 2800;
Data.Shells(8).rho = 2800;
Data.Shells(9).rho = 2800;

p  = [ 2800     ,  2800     ,  2800     ,  2800     ,  2800      , ... 
       2800     ,  2800     ,  2800     ,  2800     ]; 

% Static friction coifficient   
Data.Shells(1).mu = 0.61;
Data.Shells(2).mu = 0.61;
Data.Shells(3).mu = 0.61;
Data.Shells(4).mu = 0.61;
Data.Shells(5).mu = 0.61;
Data.Shells(6).mu = 0.61;
Data.Shells(7).mu = 0.61;
Data.Shells(8).mu = 0.61;
Data.Shells(9).mu = 0;

uu = [ 0.61     ,  0.61     ,  0.61     ,  0.61     ,  0.61      , ... 
       0.61     ,  0.61     ,  0.61     ,  0        ]; 

% Shear modulus  
Data.Shells(1).G.rz = 28*10^9;
Data.Shells(2).G.rz = 28*10^9;
Data.Shells(3).G.rz = 28*10^9;
Data.Shells(4).G.rz = 28*10^9;
Data.Shells(5).G.rz = 28*10^9;
Data.Shells(6).G.rz = 28*10^9;
Data.Shells(7).G.rz = 28*10^9;
Data.Shells(8).G.rz = 28*10^9;
Data.Shells(9).G.rz = 28*10^9;

G_rz = [ 28*10^9  ,  28*10^9  ,  28*10^9  ,  28*10^9  ,  28*10^9   , ...
       28*10^9  ,  28*10^9  ,  28*10^9  ,  28*10^9        ];      
% Cost  
Data.Shells(1).Cost = 50;
Data.Shells(2).Cost = 50;
Data.Shells(3).Cost = 50;
Data.Shells(4).Cost = 50;
Data.Shells(5).Cost = 50;
Data.Shells(6).Cost = 50;
Data.Shells(7).Cost = 50;
Data.Shells(8).Cost = 50;
Data.Shells(9).Cost = 50;

C = [ 50  ,  50  ,  50  ,  50  ,  50   , ...
       50  ,  50  ,  50  ,  50        ];   

% Additional input paramaters
Data.height = 0.3;
Data.rot_speed.max = 8000;
Data.rot_speed.min = 3000;
Data.stepsize = 0.000001;

hh = 0.3 ;                              % rotor height
n = 8000;                              % [rpm] rotational speed
n_ = 3000;                             % Minspeed
step  =  0.000001   ;  % [m] Radial step size used in calculations 
%Must be same as in Main_calc

% Preallocating vectors
Data.nShells = length(Data.Shells);
Data.nRadial = numel(Data.Shells(1).d.i/r_scale:Data.stepsize:Data.Shells(Data.nShells).d.o/r_scale);

m = length(ri);
r = min(ri):step:max(ro); % radial vector

% Fill idat with input data for easy transport into subroutines
% scaled to mm, diameters and MPa
idat(1,1:1:m) = ri*2000;
idat(2,1:1:m) = ro*2000;
idat(3,1:1:m) = Ec/10^9;
idat(4,1:1:m) = Er/10^9;
idat(5,1:1:m) = v;
idat(6,1:1:m) = p;
idat(7,1:1:m) = uu;
idat(9,1:1:m) = G_rz/10^9;
idat(10,1:1:m) = C;
idat(8,1)     = length(r);
idat(8,2)     = hh*1000;
idat(8,3)     = n;
idat(8,4)     = n_;
idat(8,5)     = 9;

end

