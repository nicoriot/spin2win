function [ idat ] = import_n_format( impraw, begline )
% Syntax:
% [Formated output data] =
% import_n_format(Imported raw data, Begining of run line markers)
 
% Format raw data into idat format for the GUI tables

% Read first seven lines after run marker
% Conataning geometrical data
for i=1:1:9
k=1;
[curr, rest]=strtok(impraw(begline+1+i));
[curr, rest]=strtok(rest);       
   while ~isempty(cell2mat(curr))
     idat(i,k) = str2double(cell2mat(curr));
     k=k+1;
     [curr, rest]=strtok(rest);
   end
end

% Cost fix
idat(10,:)=idat(9,:);

% Shear modulus fix
idat(9,:)=idat(8,:);

% clear basic settings vector
idat(8,:)=0;

% Read next five lines contaning settings data
for i=1:1:5
   [~ , rest]=strtok(impraw(begline+10+i));
   [curr]=strtok(rest);
   idat(8,i)=str2double(cell2mat(curr));
end

% Scaling to mm, diameters and MPa
idat(1,:)=idat(1,:)*2000;
idat(2,:)=idat(2,:)*2000;

idat(3,:)=idat(3,:)/10^9;
idat(4,:)=idat(4,:)/10^9;
idat(9,:)=idat(9,:)/10^9;

idat(8,2)=idat(8,2)*1000;
end

