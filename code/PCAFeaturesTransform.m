function [newTraFea, newTesFea] = PCAFeaturesTransform(traFea,tesFea,~,thresh)
if ~nargin, unitTest; return; end
if ~exist('thresh','var') || isempty(thresh), thresh = 95; end
%compute principal components and weights
[D,~,~,~,explained] = pca(traFea);
X = traFea*D';
%find a number of pc that explain 95% of the variance in the data
thrIdx = find(cumsum(explained)>=thresh,1);
%transformed training features
newTraFea = X(:,1:thrIdx)*D(1:thrIdx,1:thrIdx);
%transform test features using the dictionary learned on the
%training features
XTes = tesFea*D';                   %pca coefficients
newTesFea = XTes(:,1:thrIdx)*D(1:thrIdx,1:thrIdx);

function unitTest
clear, clc, close all
load fisheriris
nDim = 4;
nObs = 150;
traFea = meas(1:nObs,1:nDim);
tesFea = traFea;
traCat = species(1:nObs);
traFea = traFea-repmat(mean(traFea),nObs,1);
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
xlabel('Sepal length');
ylabel('Sepal width');