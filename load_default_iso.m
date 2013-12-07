function [Data] = load_default_iso()
% Returns a set of data that can be used as input data in Spin2Win

% Inner diameter [mm]
r_scale = 2000; % converts from radius in [m] to diameter in [mm]
Data.Shells(1).d.i = 0.100*r_scale;
Data.Shells(2).d.i = 0.125*r_scale;
Data.Shells(3).d.i = 0.150*r_scale;
Data.Shells(4).d.i = 0.175*r_scale;
Data.Shells(5).d.i = 0.200*r_scale;
Data.Shells(6).d.i = 0.225*r_scale;
Data.Shells(7).d.i = 0.250*r_scale;
Data.Shells(8).d.i = 0.275*r_scale;
Data.Shells(9).d.i = 0.300*r_scale;

% Outer diameter [mm]
Data.Shells(1).d.o = 0.12525*r_scale;
Data.Shells(2).d.o = 0.15025*r_scale;
Data.Shells(3).d.o = 0.17525*r_scale;
Data.Shells(4).d.o = 0.20025*r_scale;
Data.Shells(5).d.o = 0.22525*r_scale;
Data.Shells(6).d.o = 0.25025*r_scale;
Data.Shells(7).d.o = 0.27525*r_scale;
Data.Shells(8).d.o = 0.30025*r_scale;
Data.Shells(9).d.o = 0.32500*r_scale;

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

% Possions ratio, cirumference - radial direction   
Data.Shells(1).v = 0.33;
Data.Shells(2).v = 0.33;
Data.Shells(3).v = 0.33;
Data.Shells(4).v = 0.33;
Data.Shells(5).v = 0.33;
Data.Shells(6).v = 0.33;
Data.Shells(7).v = 0.33;
Data.Shells(8).v = 0.33;
Data.Shells(9).v = 0.33;

% Material density [kg/m3]
Data.Shells(1).rho = 2800;
Data.Shells(2).rho = 2800;
Data.Shells(3).rho = 2800;
Data.Shells(4).rho = 2800;
Data.Shells(5).rho = 2800;
Data.Shells(6).rho = 2800;
Data.Shells(7).rho = 2800;
Data.Shells(8).rho = 2800;
Data.Shells(9).rho = 2800;

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

% Shear modulus [MPa]
Data.Shells(1).G.rz = 28;
Data.Shells(2).G.rz = 28;
Data.Shells(3).G.rz = 28;
Data.Shells(4).G.rz = 28;
Data.Shells(5).G.rz = 28;
Data.Shells(6).G.rz = 28;
Data.Shells(7).G.rz = 28;
Data.Shells(8).G.rz = 28;
Data.Shells(9).G.rz = 28;

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

% Additional input paramaters
Data.height = 0.3*1e3; %[mm]
Data.rot_speed.max = 8000; %[rpm]
Data.rot_speed.min = 3000; %[rpm]
Data.stepsize = 0.000001; 
Data.nShells = length(Data.Shells);
Data.nRadial = numel(Data.Shells(1).d.i/r_scale:Data.stepsize:Data.Shells(Data.nShells).d.o/r_scale);
Data.r = (Data.Shells(1).d.i/r_scale):Data.stepsize:(Data.Shells(end).d.o/r_scale);
Data.figurenr = 1;
end

