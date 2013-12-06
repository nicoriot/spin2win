% loading and plotting strains from mout.mat
% dummy file
% to be expanded and integrated into plotter.m at later date,
% when strain calcualtions work fully.

load('mout.mat')
figure(4); 
plot(mout.r*2000,mout.e_r(2,:)) 
hold on
zz(1:length(mout.r))=mout.e_z(2,1);
plot(mout.r*2000,zz)
plot(mout.r*2000,mout.e_c(2,:))