function [ idat ] = load_default_ort()
% Syntax:
% [Output data] = load_default_ort()

% Generates a set of data that can be used as input data in Spin2Win

% Input paramaters
%                  CARBON-COMPOSITE SETTINGS
%       C1          C2          C3          C4         C5 
%       C6          C7          C8          C9
% [m] Inner radius
ri = [ 0.100    ,  0.125    ,  0.150    ,  0.175    ,  0.200     , ... 
       0.225    ,  0.250    ,  0.275    ,  0.300    ]; 
% [m] Outer radius   
ro = [ 0.12525  ,  0.15025  ,  0.17525  ,  0.20025  ,  0.22525   , ...
       0.25025  ,  0.27525  ,  0.30025  ,  0.32500  ]; 
% [Pa] Youngs Modulus, Hoop   
Ec = [ 155*10^9 ,  155*10^9 ,  155*10^9 ,  155*10^9 ,  155*10^9  , ...
       155*10^9 ,  155*10^9 ,  155*10^9 ,  155*10^9 ]; 
% [Pa ]Youngs Modulus, Radial   
Er = [ 9*10^9   ,  9*10^9   ,  9*10^9   ,  9*10^9   ,  9*10^9    , ...
       9*10^9   ,  9*10^9   ,  9*10^9   ,  9*10^9   ]; 
% Possions ratio, cirumference - radial direction   
v  = [ 0.25     ,  0.25     ,  0.25     ,  0.25     ,  0.25      , ...
       0.25     ,  0.25     ,  0.25     ,  0.25     ]; 
% [kg/m3] Material density   
p  = [ 1600     ,  1600     ,  1600     ,  1600     ,  1600      , ...
       1600     ,  1600     ,  1600     ,  1600     ]; 
% Static friction coifficient   
uu = [ 0.61     ,  0.61     ,  0.61     ,  0.61     ,  0.61      , ...
       0.61     ,  0.61     ,  0.61     ,  0        ]; 
% Shear modulus  
G_rz = [ 3.6*10^9 , 3.6*10^9 , 3.6*10^9 , 3.6*10^9 , 3.6*10^9 , ...
         3.6*10^9 , 3.6*10^9 , 3.6*10^9 , 3.6*10^9  ];    

% Adittional input paramaters
hh = 0.3 ;                              % rotor height
n = 10000;                              % [rpm] rotational speed
n_ = 5000;                              % Minimal speed

% Settings
step  =  0.000001   ;  % [m] Radial step size used in calculations 
% Must be same as in Main_calc

% Preallocating vectors, mostly diffrent constants used in loops
m = length(ri);
r = min(ri):step:max(ro); % radial vector

% Fill idat with input data for easy transport into subroutines
% Scaled to mm, diameters and MPa
idat(1,1:1:m) = ri*2000;
idat(2,1:1:m) = ro*2000;
idat(3,1:1:m) = Ec/10^9;
idat(4,1:1:m) = Er/10^9;
idat(5,1:1:m) = v;
idat(6,1:1:m) = p;
idat(7,1:1:m) = uu;
idat(9,1:1:m) = G_rz/10^9;
idat(8,1)     = length(r);
idat(8,2)     = hh*1000;
idat(8,3)     = n;
idat(8,4)     = n_;
idat(8,5)     = 9;

end

