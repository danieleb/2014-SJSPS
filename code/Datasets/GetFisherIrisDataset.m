function [fea,cat] = GetFisherIrisDataset(visu,dim)
% Gets the fisher iris dataset from the web and returns the following variables
% - fea: matrix containing sepal length, width, petal length and width along the columns
% - cat: cell containing the classes
if ~exist('visu','var'), visu = false; end
url = 'http://archive.ics.uci.edu/ml/machine-learning-databases/iris/iris.data';
s = urlread(url);
C = textscan(s,'%f,%f,%f,%f,%s');
fea = [C{1} C{2} C{3} C{4}];
cat = C{5};

if visu
    if ~exist('dim','var'), dim = 3; end
    fea = fea(:,1:dim);
    figure
    if dim==2
        gscatter(fea(:,1),fea(:,2),cat);
    elseif dim==3
        gscatter3(fea(:,1),fea(:,2),fea(:,3),cat);
    end
    xlabel('sepal length')
    ylabel('sepal width')
    zlabel('petal length')
    legend('setosa','versicolor','virginica')
    [a,b] = fileparts(which(mfilename));
    savegraph([a,filesep,b]);
end

