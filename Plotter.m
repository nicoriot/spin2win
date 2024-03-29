function Plotter(Ten_c, Ten_r, intf, r_mark, r_mark2, r, ri, fignum, stamp, prefix)
% Syntax:
% Plotter( Hoop stress, Radial stress, Radial displacements, ...
% Shell markers upper half, Shell markers lower half, Radial vector, ...
% Inner radius, Figure number, Time stamp, Run prefix)

% Plot input data in one figure with three subplots.

if ishandle(fignum) == 1
close(figure(fignum))
end

%%%%%%%%%%%%%%%%%%%%%
% KURT plotter haxx %
do_it = 0;
filen = 'Kurt_L2_raw';

if do_it == 1;
    len = length(Ten_c(3,:));
    [tmp_hoop, tmp_rad] = Kurt_format(len,filen);
    Ten_c(3,:)=Ten_c(3,:)+tmp_hoop;
    Ten_r(3,:)=Ten_r(3,:)+tmp_rad;
    Ten_c(1,:)=Ten_c(1,:)+tmp_hoop;
    Ten_r(1,:)=Ten_r(1,:)+tmp_rad;
    
end
%%%%%%%%%%%%%%%%%%%%%
% Settings
step  =  0.000001   ;  % [m] Radial step size used in calculations
% Must be same as in Main_calc

% Make figure
figure(fignum)

% Set figure name to same as stamp
set(figure(fignum),'Name',[stamp, ', ', prefix])

% Plotting hoop stresses
subplot(3,1,1)
hold on
plot(r*1000,Ten_c(1,:),'--r');
plot((r_mark*step+ri(1))*1000,Ten_c(1,r_mark),'*r');
plot(((r_mark2-1)*step+ri(1))*1000,Ten_c(1,(r_mark2)),'*r');

plot(r*1000,Ten_c(2,:),'--b');
plot((r_mark*step+ri(1))*1000,Ten_c(2,r_mark),'*b');
plot(((r_mark2-1)*step+ri(1))*1000,Ten_c(2,(r_mark2)),'*b');

plot(r*1000,Ten_c(3,:),'-k');
plot((r_mark*step+ri(1))*1000,Ten_c(3,r_mark),'*k');
plot(((r_mark2-1)*step+ri(1))*1000,Ten_c(3,(r_mark2)),'*k');

% Give axis names
hold off
title 'Hoop stress [Blue: Pure centrifugal , Red: Press-fit at standstill , Black: Press-fit while rotating]'
ylabel 'Stress [MPa]'

% Plotting radial stresses
subplot(3,1,2)
hold on

plot(r*1000,Ten_r(1,:),'--r');
plot((r_mark*step+ri(1))*1000,Ten_r(1,r_mark),'*r');

plot(r*1000,Ten_r(2,:),'--b');
plot(((r_mark)*step+ri(1))*1000,Ten_r(2,r_mark),'*b');

plot(r*1000,Ten_r(3,:),'-k');
plot((r_mark*step+ri(1))*1000,Ten_r(3,r_mark),'*k');

% Give axis names
hold off
title 'Radial stress  [Blue: Pure centrifugal , Red: Press-fit at standstill , Black: Press-fit while rotating]'
ylabel 'Stress [MPa]'

% Plotting interface distrubution
subplot(3,1,3)
hold on

plot(r*1000,intf(2,:)*10^6,'b');
plot((r_mark*step+ri(1))*1000,intf(2,r_mark)*10^6,'*b');
plot(((r_mark2-1)*step+ri(1))*1000,intf(2,(r_mark2))*10^6,'*b');

plot(r*1000,intf(1,:)*10^6,'r');
plot((r_mark*step+ri(1))*1000,intf(1,r_mark)*10^6,'*r');
plot(((r_mark2-1)*step+ri(1))*1000,intf(1,(r_mark2))*10^6,'*r');

% Give axis names
hold off
title 'Radial displacement [Red: Displacement at standstill , Blue: Displacement while rotating]'
ylabel 'Displacement [�m]'
xlabel 'Radius [mm]'

end

