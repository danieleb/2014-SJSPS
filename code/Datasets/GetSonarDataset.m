function [fea,cat] = GetSonarDataset()
% Gets the sonar dataset from the web
url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/undocumented/connectionist-bench/sonar/sonar.all-data';
s = urlread(url);
format = [];
for i=1:60, format = [format, sprintf('%%f ')]; end
format = [format, '%s'];
C = textscan(s,format,'Delimiter',',');
fea = [C{1:60}];
cat = C{61};
