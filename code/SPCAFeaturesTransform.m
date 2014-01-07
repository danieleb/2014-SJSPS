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

%find a number of pc that explain thresh% of the variance in the data or
%min 2d
thrIdx = 1:find(thresh-cumsum(explained)>=100*eps); %use this to avoid machine errors
if length(thrIdx) < 2, thrIdx = [1,2]; end

%transformed training features
U = D(:,thrIdx);
newTraFea = X'*U; %return n observations along the rows of matrix
newTesFea = tesFea*U;

function newCat = NumericEquivalent(cat)
%takes a cell of string categories and returns a vector of equivalent
%integer values. Example: {'a','b','c','b','a'} => [1,2,3,2,1].
cat = cellstr(cat); %transform to string in case of nominal categories
nObs = length(cat);
uniCat = unique(cat);
nCat = length(uniCat);
newCat = zeros(nObs,1);
for i=1:nCat
    newCat(strcmp(cellstr(cat),uniCat(i))) = i-1;
end

function unitTest
clear, clc
traObsPer = 0.7;
load fisheriris
nDim = 4;
nObs = 150;
iObs = randperm(nObs);
traFea = meas(iObs(1:floor(traObsPer*nObs)),1:nDim);
tesFea = meas(iObs(floor(traObsPer*nObs)+1:end),1:nDim);
traCat = species(iObs(1:floor(traObsPer*nObs)));
tesCat = species(iObs(floor(traObsPer*nObs)+1:end));

% Center features
mFea = mean(traFea);
traFea = traFea-repmat(mFea,floor(traObsPer*nObs),1);
tesFea = tesFea-repmat(mFea,nObs-floor(traObsPer*nObs),1);

figure, subplot(1,2,1)
title('original features');
Plot2dTraTes(traFea,tesFea,traCat,tesCat);
xlabel('Sepal length');
ylabel('Sepal width');

[newTraFea,newTesFea] = SPCAFeaturesTransform(traFea,tesFea,traCat);
subplot(1,2,2)
title('Transformed features');
if size(newTraFea,2)>1
    Plot2dTraTes(newTraFea,newTesFea,traCat,tesCat);
end
xlabel('1st principal component');
ylabel('2nd principal component');


function Plot2dTraTes(traFea,tesFea,traCat,tesCat)
cg = 'rgb'; % color scheme
[sTraCat,idx] = sort(traCat);  %sort categories
gscatter(traFea(idx,1), traFea(idx,2), sTraCat,cg,'osd');
hold on
[sTesCat,idx] = sort(tesCat);
h = gscatter(tesFea(idx,1),tesFea(idx,2),sTesCat,cg,'osd');
for i=1:length(h)
    set(h(i),'MarkerFaceColor',cg(i))
end
