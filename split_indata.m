function [ geomdat, setdat ] = split_indata( senddat )
% Syntax:
% [Geometrical data, Settings data] =
% split_indata(Input raw data)

% Pre allocate variables
geomdat = zeros(9,length(senddat(1,:)));
setdat  = zeros(1,5);

% Extract geometrical data from input
for i=1:1:7
    for k=1:1:length(senddat(1,:))
    geomdat(i,k) = senddat(i,k);
    end
end

for k=1:1:length(senddat(1,:))
geomdat(8,k) = senddat(9,k);
end

for k=1:1:length(senddat(1,:))
geomdat(9,k) = senddat(10,k);
end

% Extract settings data from input
i = 8;
for k=1:1:5
    setdat(1,k) = senddat(i,k);
end

end