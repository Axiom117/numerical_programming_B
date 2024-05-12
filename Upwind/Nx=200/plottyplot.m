dd = dir('*.csv');

fileNames = {dd.name}; 

data = cell(numel(fileNames),2);
data(:,1) = regexprep(fileNames, '.csv','');

for i = 1:numel(fileNames)    
   data{i,2} = dlmread(fileNames{i}, ',', 1, 0);
end

fig=figure();
hold on;

for j = 1:numel(fileNames)

    XY = data{j,2};
    X = XY(:,1);
    Y = XY(:,2);

    
    if(j <= 7)
        plot(X,Y, 'LineWidth', 2);
    else
        plot(X,Y, '--', 'LineWidth', 2);
    end
end

legend(fileNames);
title("Differential wave Equation Solution Upwind scheme, Nx=200");