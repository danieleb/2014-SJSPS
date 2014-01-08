function [newTraFea, newTesFea] = PCAFeaturesTransform(traFea,tesFea,~,par)
if ~nargin, unitTest; return; end
if ~exist('par','var') || isempty(par), par=struct; end
def.thresh = 95;
def.subSpaRan = [];
par = setdefaultoptions(par,def);

%compute principal components and weights
[D,~,~,~,explained] = pca(traFea);
if par.subSpaRan
    D = D(:,1:par.subSpaRan);
else
%find a number of pc that explain 95% of the variance in the data
    thrIdx = find(cumsum(explained)>=par.thresh,1);
    D = D(:,1:thrIdx);
end
Project = @(Phi,x) Phi*pinv(Phi)*x;
%transformed training features
newTraFea = Project(D,traFea);
newTesFea = Project(D,tesFea);

function unitTest
clear, clc
load fisheriris
nDim = 4;
nObs = 150;
traFea = meas(1:nObs,1:nDim);
tesFea = traFea;
traCat = species(1:nObs);
traFea = traFea-repmat(mean(traFea),nObs,1);
figure
gscatter(traFea(:,1), traFea(:,2), traCat,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');

[newTraFea,newTesFea] = PCAFeaturesTransform(traFea,tesFea);
size(newTraFea,2);
figure
if size(newTraFea,2)<2
    gscatter(newTraFea(:,1), newTraFea(:,1), traCat,'rgb','osd');
else
    gscatter(newTraFea(:,1), newTraFea(:,2), traCat,'rgb','osd');
end
xlabel('1st principal component');
ylabel('2nd principal component');