function [ senddat ] = merge_indata( geomdat, setdat )
% Syntax:
% [Output matrix] =
% merge_indata(Geometrical data, Settings data)

% Merges geometrical and settings data into a singel matrix
% for easy transport between subfunctions

% Pre allocate matrix
senddat = zeros(10,length(geomdat(1,:)));

% Copy geometrical data into output matrix
for i=1:1:7
    for k=1:1:length(geomdat(1,:))
    senddat(i,k) = geomdat(i,k);
    end
end

for k=1:1:length(geomdat(1,:))
senddat(9,k) = geomdat(8,k);
end

for k=1:1:length(geomdat(1,:))
senddat(10,k) = geomdat(9,k);
end


% Copy settings data into output matrix
i = 8;
for k=2:1:4
    senddat(i,k) = setdat(k,1);
end

% Insert number of active shells into output matrix
senddat(i,5) = setdat(1,1);

end


