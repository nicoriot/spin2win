function [texten, Dat_line] = ReadFile(infile)
% Syntax:
% [Loaded text, Beginning of run line marks] =
% ReadFile(File name of old log file)

% Check if file name was passed on correctly
% Return -1 if load failed
if infile == 0
    texten = -1;
    Dat_line = -1;
else

% Open and read file
fid = fopen(infile);
k = 1;

% Scan each line from file and save
while ~feof(fid)
    tline = fgetl(fid);
    texten(k,1) = cellstr(tline);
    k=k+1;   
end

% Close file
fclose(fid);

% Scan for 'Date and time of run:' string and mark each hit with
% corresponding line number
al = 'Date and time of run:';
i = 1;
for k=1:1:length(texten)
    A = sscanf(char(texten(k)), '%c', [1 21]);
    if strcmp(A, al)
        Dat_line(i) = k;
        i =i+1;
    end
end

%Swap Dat_line direction
Dat_line(length(Dat_line):-1:1) = Dat_line;
end

end



        