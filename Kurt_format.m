function [ WW_hoop, WW_rad] = Kurt_format(leng, filename)
% Formats WetWind, SwereaSicomp output data to correct length
% for use in spin 2 win plots

dat=load(filename);

dat= struct2array(dat(1));

r = min(dat(:,1))/1000:0.000001:max(dat(:,1)/1000);

rad = interp1(dat(:,1)/1000,dat(:,2),r);
hoop = interp1(dat(:,1)/1000,dat(:,3),r);

diff = leng - length(rad) + 1;

WW_hoop = zeros(1,leng);
WW_rad  = zeros(1,leng);

WW_hoop(diff:leng) = WW_hoop(diff:leng)+hoop(1:length(hoop));
WW_rad(diff:leng) = WW_rad(diff:leng)+rad(1:length(rad));

end

