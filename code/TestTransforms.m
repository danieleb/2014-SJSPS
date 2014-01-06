function res = TestTransforms(traFea,tesFea,traCat,tesCat,par)

addpath(genpath(pwd));
reset(RandStream.getGlobalStream)

%% Parameters and defaults
if ~exist('par','var') || isempty(par), par = struct; end

def.nAto = 10;
def.subSpaRan = 1;
def.perActAto = 50;
def.nKNNs = 5;

par = setdefaultoptions(par,def);

if ~exist('traFea','var') || isempty(traFea) || ~exist('tesFea','var') || isempty(tesFea)...
        || ~exist('traCat','var') || isempty(traCat) || ~exist('tesCat','var') || isempty(tesCat)
    [meas,species] = GetFisherIrisDataset;
    nDim = 3;
    nObs = 150;
    v = randperm(nObs);
    traFea = meas(v(1:100),1:nDim);
    tesFea = meas(v(101:150),1:nDim);
    traCat = species(v(1:100));
    tesCat = species(v(101:150));
    feaMea = mean(traFea);
    feaStd = std(traFea);
    traFea = (traFea-repmat(feaMea,100,1))./repmat(feaStd,100,1);
    tesFea = (tesFea-repmat(feaMea,50,1))./repmat(feaStd,50,1);
end

%% Classify using linear discriminant analysis
cvPar = cvpartition([traCat;tesCat],'kfold',5);                 %create 5fold partition

conMatFun = @(xTra,yTra,xTes,yTes) (confusionmat(yTes,...               %Compute confusion matrix
    nominal(knnclassify(xTes,xTra,yTra,par.nKNNs))));
conMat = crossval(conMatFun,[traFea;tesFea],nominal([traCat;tesCat]),'partition',cvPar);

nCat = length(unique(nominal([traCat;tesCat])));
conMat = reshape(sum(conMat),nCat,nCat);
mcr = sum(sum(conMat-diag(diag(conMat))))/sum(sum(conMat));

%% Test PCA feature transform
[pcaTraFea, pcaTesFea] = PCAFeaturesTransform(traFea,tesFea);

pcaConMat = crossval(conMatFun,[pcaTraFea;pcaTesFea],nominal([traCat;tesCat]),'partition',cvPar);
pcaConMat = reshape(sum(pcaConMat),nCat,nCat);
pcaMcr = sum(sum(pcaConMat-diag(diag(pcaConMat))))/sum(sum(pcaConMat));

%% Test LDA feature transform
[ldaTraFea,ldaTesFea] = LDAFeaturesTransform(traFea,tesFea,traCat);

conMatFun = @(xTra,yTra,xTes,yTes) (confusionmat(yTes,...               %Compute confusion matrix
    nominal(knnclassify(xTes,xTra,yTra,par.nKNNs))));
ldaConMat = crossval(conMatFun,[ldaTraFea;ldaTesFea],nominal([traCat;tesCat]),'partition',cvPar);
ldaConMat = reshape(sum(ldaConMat),nCat,nCat);
ldaMcr = sum(sum(ldaConMat-diag(diag(ldaConMat))))/sum(sum(ldaConMat));

%% Test IPR feature transform
[~, iprTraFea, iprTesFea] = IPRLearnDisSub(traFea,tesFea,traCat,par);

% conMatFun = @(xTra,yTra,xTes,yTes) (confusionmat(yTes,...               %Compute confusion matrix
%     nominal(IPRClassify(xTes,xTra,yTra,traSubSpa,uniCat))));
conMatFun = @(xTra,yTra,xTes,yTes) (confusionmat(yTes,...
    nominal(knnclassify(xTes,xTra,yTra,par.nKNNs))));
iprConMat = crossval(conMatFun,[iprTraFea;iprTesFea],nominal([traCat;tesCat]),'partition',cvPar);
iprConMat = reshape(sum(iprConMat),nCat,nCat);
iprMcr = sum(sum(iprConMat-diag(diag(iprConMat))))/sum(sum(iprConMat));


%% Plots
J1 = @(fea,cat) trace(pinv(ClassificationDiscriminant.fit(fea,cat).Sigma)*ClassificationDiscriminant.fit(fea,cat).BetweenSigma);
J2 = @(fea,cat) trace(ClassificationDiscriminant.fit(fea,cat).BetweenSigma)/trace(ClassificationDiscriminant.fit(fea,cat).Sigma);
J3 = @(fea,cat) log(abs(det(ClassificationDiscriminant.fit(fea,cat).BetweenSigma))) - ClassificationDiscriminant.fit(fea,cat).LogDetSigma;

mets = {'none','PCA','LDA','IPR'};
nMet = length(mets);
sqrNMet = ceil(sqrt(nMet));
% Display features
close all
for iMet=1:nMet
    switch iMet
        case 1
            fea = traFea;
            mat = conMat;
            err = mcr;
        case 2
            fea = [pcaTraFea zeros(size(pcaTraFea,1),1)];
            mat = pcaConMat;
            err = pcaMcr;
        case 3
            fea = ldaTraFea;
            mat = ldaConMat;
            err = ldaMcr;
        case 4
            fea = iprTraFea;
            mat = iprConMat;
            err = iprMcr;
    end
   figure(1), subplot(sqrNMet,sqrNMet,iMet)
   ScatterFeatures(fea,traCat,mets{iMet});
%    figure(2), subplot(sqrNMet,sqrNMet,iMet)
%    DispConfMat(mat,unique([traCat;tesCat])), colormap jet;
%    title([mets{iMet} ' features transforms. MCR: ' num2str(err)])
end

res.mcr = [mcr pcaMcr ldaMcr iprMcr]';
res.mets = mets;
res.traFea = {traFea,pcaTraFea,ldaTraFea,iprTraFea};
res.tesFea = {tesFea,pcaTesFea,ldaTesFea,iprTesFea};
end

function ScatterFeatures(traFea,traCat,method)
traFea = NormalizeFeatures(traFea);
if size(traFea,2)<2
    gscatter(traFea(:,1), traFea(:,1), traCat,'rgb','osd');
else
    uniCat = unique(traCat);
    col = {'r','g','b','y','k','c'};
    mar = {'o','s','d','x','^','>'};
    for iCat=1:length(uniCat)
        ind = strcmp(uniCat(iCat),traCat);
        scatter3(traFea(ind,1),traFea(ind,2),traFea(ind,3),100,col{iCat},'Marker',mar{iCat});
        hold on
    end
    
end
title(['features transform: ' method]);
xlabel('Sepal length')
ylabel('Sepal width')

end