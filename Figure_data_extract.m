% Figure data extractor for subplotted figures
% custom made to extract Spin2Win figure data.

% Currently standalone, to be implemented into a data extraction tool.

% get open figure handle number
hand = gcf;

% Get all data
alldat = findall(hand,'Type','line');

% Get all x-data
x_all=get(alldat,'Xdata');

% Get all y-data
y_all=get(alldat,'YData');

% Find double locations
i = 1;
for k=[1:length(x_all)]
tmp = x_all(k,:);
tmp = cell2mat(tmp);
if length(tmp)>100
    x_spot(i)=k;
    i=i+1;
end
end

i = 1;
y_spot = 0;
for k=[1:length(y_all)]
tmp = y_all(k,:);
tmp = cell2mat(tmp);
if length(tmp)>100
    y_spot(i)=k;
    i=i+1;
end
end

% Extract plot data
for k=[1:length(x_spot)]
x_cell(1,k) = (x_all(x_spot(k),:));
x_dat(k,:)=cell2mat(x_cell(k));
end

for k=[1:length(y_spot)]
y_cell(1,k) = (y_all(y_spot(k),:));
y_dat(k,:)=cell2mat(y_cell(k));
end

