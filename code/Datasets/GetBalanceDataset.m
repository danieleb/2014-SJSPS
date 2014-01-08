function [fea,cat] = GetBalanceDataset()
% Gets the fisher balance dataset from the web (see http://archive.ics.uci.edu/ml/datasets/Balance+Scale)
url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/balance-scale/balance-scale.data';
s = urlread(url);
C = textscan(s,'%s %f %f %f %f','Delimiter',',');
fea = [C{2:5}];
cat = C{1};
