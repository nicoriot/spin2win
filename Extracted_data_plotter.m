% Plot imported figure data data
% given line numbers, imported data
% is counting backwards.

% Rename data 1 = SW data, 2 = S2W data

% Hoop extractor
Ten_c1 = Ydat_SW(3,:);
r1     = Xdat_SW(3,:);

Ten_c2 = Ydat_S2W(3,:);
r2     = Xdat_S2W(3,:);

% Radial extractor
Ten_r1 = Ydat_SW(2,:);
r3     = Xdat_SW(2,:);

Ten_r2 = Ydat_S2W(2,:);
r4     = Xdat_S2W(2,:);

% R-disp extractor
intf1 = Ydat_SW(1,:);
r5     = Xdat_SW(1,:);

intf2 = Ydat_S2W(1,:);
r6     = Xdat_S2W(1,:);

% Plotting hoop stresses
subplot(3,1,1)
hold on
plot(r1,Ten_c1,'--r');

plot(r2,Ten_c2,'b');

% Give axis names
hold off
legend('Layup 2 standard','Layup 2 50-50')
title 'Hoop stress'
ylabel 'Stress [MPa]'

% Plotting radial stresses
subplot(3,1,2)
hold on

plot(r3,Ten_r1,'--r');

plot(r4,Ten_r2,'b');

% Give axis names
hold off
title 'Radial stress'
ylabel 'Stress [MPa]'

% Plotting interface distrubution
subplot(3,1,3)
hold on

plot(r5,intf1,'--r');

plot(r6,intf2,'b');

% Give axis names
hold off
title 'Radial displacement'
ylabel 'Displacement [�m]'
xlabel 'Radius [mm]'