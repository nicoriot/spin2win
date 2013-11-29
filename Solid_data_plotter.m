% Solid data plotter and axis fixer
% All data must be extracted forom the soliworks model using the same
% probe sensor

dh = Rotor0x2D30k5050_hoop   ; % hoop data matrix name
dr = Rotor0x2D30k5050_radial ; % radial data matrix name
di = Rotor0x2D30k5050_rdisp  ; % radial displacemnt matrix name
n = 1001                   ;

% sign fix
dh(:,3)=abs(dh(:,3));
dr(:,3)=abs(dr(:,3));
di(:,3)=abs(di(:,3));

% ri ro fix
ri = round(min(dh(:,3)));
ro = round(max(dh(:,3)));

% dm fix
ld = length(dh(:,1)) ;
dh(1,3) = ri ;
dh(ld,3) = ro ;

%extract data x axis
dx = dh(:,3) ;

%extract data y axis, hoop
dyh = dh(:,2) ;

%extract data y axis, radial
dyr = dr(:,2) ;

%extract data y axis, disp
dyi = di(:,2) ;

% make x axis
x = linspace(ri,ro,n) ;

% interpolate y axis, hoop
yh = interp1(dx,dyh,x)  ;

% interpolate y axis, radial
yr = interp1(dx,dyr,x)  ;

% interpolate y axis, radial
yi = interp1(dx,dyi,x)  ;

% Plot hoop
figure(5)
subplot(3,1,1)
plot(x,yh)
ylabel('Stress [MPa]')
title('Hoop stress')
xlim([ri-2 ro+2])
ylim([-50 650])

% Plot radial
subplot(3,1,2)
plot(x,yr)
ylabel('Stress [MPa]')
title('Radial stress')
xlim([ri-2 ro+2])
ylim([-50 10])

% Plot rdisp
subplot(3,1,3)
plot(x,yi)
xlabel('Radius [mm]')
ylabel('Displacement [µm]')
title('Radial displacement')
xlim([ri-2 ro+2])
ylim([0 900])