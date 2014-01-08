function [fea,cat] = GetFisherIrisDataset()
% Gets the fisher iris dataset from the web and returns the following variables
% - fea: matrix containing sepal length, width, petal length and width along the columns
% - cat: cell containing the classes
url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data';
s = urlread(url);
C = textscan(s,'%f,%f,%f,%f,%s');
fea = [C{1} C{2} C{3} C{4}];
cat = C{5};