function [fea,cat] = GetParkinsonsDataset()
% Gets the Parkinsons dataset from the web
url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/parkinsons/parkinsons.data';
s = urlread(url);
format = '%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %s %f %f %f %f %f %f';
C = textscan(s(230:end),format,'Delimiter',',');
fea = [C{2:17} C{19:24}];
cat = C{18};
