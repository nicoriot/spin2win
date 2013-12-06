function [ senddat,d ] = merge_indata( geomdat, setdat )
% Syntax:
% [Output matrix] =
% merge_indata(Geometrical data, Settings data)

n = length(geomdat(1,:));
senddat = zeros(10,n);
senddat(1:7,1:n) = geomdat(1:7,1:n);
senddat(9,1:n) = geomdat(8,1:n);
senddat(10,1:n) = geomdat(9,1:n);
senddat(8,2:4) = setdat(2:4,1);
senddat(8,5) = setdat(1,1);

d.m = 0;

end


