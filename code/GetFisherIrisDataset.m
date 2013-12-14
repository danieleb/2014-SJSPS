function [fea,cat] = GetFisherIris()
% Gets the fisher iris dataset from the web and returns the following variables
% - fea: matrix containing sepal length, width, petal length and width along the columns
% - cat: cell containing the classes
url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data';
s = urlread(url);
[seplen,sepwid,petlen,petwid,cat] = strread(s,"%f,%f,%f,%f,%s");
fea = [seplen, sepwid, petlen, petwid];
