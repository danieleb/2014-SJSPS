function [newTraFea, newTesFea] = SPCAFeaturesTransform(traFea,tesFea,cat,thresh)
%supervised pca implementation. For details, please refer to Barshan et al.
%, Pattern Recognition, 44(2011) 1357--1371

if ~nargin, unitTest; return; end
if ~exist('thresh','var') || isempty(thresh), thresh = 100; end

%compute centering matrix
X = traFea';
n = size(X,2); %number of observations
H = eye(n) - ones(n,1)*ones(n,1)'/n;

%compute objective matrix Q
if ~isnumeric(cat), cat = NumericEquivalent(cat); end %if category is not numeric, compute the equivalent
L = cat*cat'; %kernel of responses (classes)
Q = X*H*L*H*X';

%compute principal components and weights
[D,~,explained] = pcacov(Q); %use pcacov because Q is a covariance matrix

%find a number of pc that explain thresh% of the variance in the data
thrIdx = find(cumsum(explained)>=thresh,1);

%transformed training features
U = D(:,1:thrIdx);
newTraFea = X'*U; %return n observations along the rows of matrix
newTesFea = tesFea*U;

function newCat = NumericEquivalent(cat)
%takes a cell of string categories and returns a vector of equivalent
%integer values. Example: {'a','b','c','b','a'} => [1,2,3,2,1].
nObs = length(cat);
uniCat = unique(cat);
nCat = length(uniCat);
newCat = zeros(nObs,1);
for i=1:nCat
    newCat(strcmp(cat,uniCat(i))) = i-1;
end

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

newTraFea = SPCAFeaturesTransform(traFea,tesFea,traCat);
figure
if size(newTraFea,2)<2
    gscatter(newTraFea(:,1), newTraFea(:,1), traCat,'rgb','osd');
else
    gscatter(newTraFea(:,1), newTraFea(:,2), traCat,'rgb','osd');
end
xlabel('1st principal component');
ylabel('2nd principal component');