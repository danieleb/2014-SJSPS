function res = TestTransforms(traFea,tesFea,traCat,tesCat,par)

addpath(genpath(pwd)); %add path to matlab search path
reset(RandStream.getGlobalStream) %reset random stream to reproduce results exactly

%% Parameters and defaults
if ~exist('par','var') || isempty(par), par = struct; end

def.nAto = 10; %number of atoms
def.subSpaRan = 1; %subspace range
def.perActAto = 50; %percentage of active atoms
def.nKNNs = 5; %number of nearest neighbours for classification

par = setdefaultoptions(par,def);

if ~exist('traFea','var') || isempty(traFea) || ~exist('tesFea','var') || isempty(tesFea)...
        || ~exist('traCat','var') || isempty(traCat) || ~exist('tesCat','var') || isempty(tesCat)
    [meas,species] = GetFisherIrisDataset; %get data from iris dataset
    nDim = 3; %data has 4 dimensions
    nObs = 150;  %and 150 observations
    v = randperm(nObs); %select random permutation of observations
    perTraObs = 70; %percentage of observations used for training
    nTraObs = floor(nObs*perTraObs/100); %number of training observations
    traFea = meas(v(1:nTraObs),1:nDim); %training data
    tesFea = meas(v(nTraObs+1:end),1:nDim); %test data
    traCat = species(v(1:nTraObs)); %training categories
    tesCat = species(v(nTraObs+1:end)); %test categories
    feaMea = mean(traFea); %mean of training categories
    feaStd = std(traFea); %standard deviation of training categories
    traFea = (traFea-repmat(feaMea,nTraObs,1))./repmat(feaStd,nTraObs,1); %normalize training data
    tesFea = (tesFea-repmat(feaMea,nObs-nTraObs,1))./repmat(feaStd,nObs-nTraObs,1); %normalize test data
end

%% Classify using original features
cvPar = cvpartition([traCat;tesCat],'kfold',5);                 %create 5fold partition

conMatFun = @(traFea,traCat,tesFea,tesCat)TransformAndClassify(traFea,traCat,tesFea,tesCat,par);
conMat = crossval(conMatFun,[traFea;tesFea],nominal([traCat;tesCat]),'partition',cvPar); %compute crossvalidation results on original features

nCat = length(unique(nominal([traCat;tesCat]))); %number of categories
conMat = reshape(sum(conMat),nCat,nCat); %confusion matrix
mcr = sum(sum(conMat-diag(diag(conMat))))/sum(sum(conMat)); %misclassification ratio

%% Test PCA feature transform
[pcaTraFea, pcaTesFea] = PCAFeaturesTransform(traFea,tesFea); %compute PCA feature transform
conMatFun = @(traFea,traCat,tesFea,tesCat)PCATransformAndClassify(traFea,traCat,tesFea,tesCat,par);
pcaConMat = crossval(conMatFun,[pcaTraFea;pcaTesFea],nominal([traCat;tesCat]),'partition',cvPar);
pcaConMat = reshape(sum(pcaConMat),nCat,nCat);
pcaMcr = sum(sum(pcaConMat-diag(diag(pcaConMat))))/sum(sum(pcaConMat));

%% Test SPCA feature transform
[spcaTraFea, spcaTesFea] = SPCAFeaturesTransform(traFea,tesFea,traCat);
conMatFun = @(traFea,traCat,tesFea,tesCat)SPCATransformAndClassify(traFea,traCat,tesFea,tesCat,par);
spcaConMat = crossval(conMatFun,[spcaTraFea;spcaTesFea],nominal([traCat;tesCat]),'partition',cvPar);
spcaConMat = reshape(sum(spcaConMat),nCat,nCat);
spcaMcr = sum(sum(spcaConMat-diag(diag(spcaConMat))))/sum(sum(spcaConMat));

%% Test LDA feature transform
[ldaTraFea,ldaTesFea] = LDAFeaturesTransform(traFea,tesFea,traCat);
conMatFun = @(traFea,traCat,tesFea,tesCat)LDATransformAndClassify(traFea,traCat,tesFea,tesCat,par);
ldaConMat = crossval(conMatFun,[ldaTraFea;ldaTesFea],nominal([traCat;tesCat]),'partition',cvPar);
ldaConMat = reshape(sum(ldaConMat),nCat,nCat);
ldaMcr = sum(sum(ldaConMat-diag(diag(ldaConMat))))/sum(sum(ldaConMat));

%% Test IPR feature transform
[~, iprTraFea, iprTesFea] = IPRLearnDisSub(traFea,tesFea,traCat,par);
conMatFun = @(traFea,traCat,tesFea,tesCat)IPRTransformAndClassify(traFea,traCat,tesFea,tesCat,par);
iprConMat = crossval(conMatFun,[iprTraFea;iprTesFea],nominal([traCat;tesCat]),'partition',cvPar);
iprConMat = reshape(sum(iprConMat),nCat,nCat);
iprMcr = sum(sum(iprConMat-diag(diag(iprConMat))))/sum(sum(iprConMat));


%% Plots
% J1 = @(fea,cat) trace(pinv(ClassificationDiscriminant.fit(fea,cat).Sigma)*ClassificationDiscriminant.fit(fea,cat).BetweenSigma);
% J2 = @(fea,cat) trace(ClassificationDiscriminant.fit(fea,cat).BetweenSigma)/trace(ClassificationDiscriminant.fit(fea,cat).Sigma);
% J3 = @(fea,cat) log(abs(det(ClassificationDiscriminant.fit(fea,cat).BetweenSigma))) - ClassificationDiscriminant.fit(fea,cat).LogDetSigma;

mets = {'none','PCA','SPCA','LDA','IPR'};
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
            tesFea = [pcaTesFea zeros(size(pcaTesFea,1),1)];
            mat = pcaConMat;
            err = pcaMcr;
        case 3
            fea = [spcaTraFea zeros(size(spcaTraFea,1),1)];
            tesFea = [spcaTesFea zeros(size(spcaTesFea,1),1)];
            mat = spcaConMat;
            err = spcaMcr;
        case 4
            fea = ldaTraFea;
            tesFea = ldaTesFea;
            mat = ldaConMat;
            err = ldaMcr;
        case 5
            fea = iprTraFea;
            tesFea = iprTesFea;
            mat = iprConMat;
            err = iprMcr;
    end
   figure(1), subplot(sqrNMet,sqrNMet,iMet)
   ScatterFeatures(fea,tesFea,traCat,tesCat,mets{iMet});
%    figure(2), subplot(sqrNMet,sqrNMet,iMet)
%    DispConfMat(mat,unique([traCat;tesCat])), colormap jet;
%    title([mets{iMet} ' features transforms. MCR: ' num2str(err)])
end

res.mcr = [mcr pcaMcr spcaMcr ldaMcr iprMcr]';
res.mets = mets;
res.traFea = {traFea,pcaTraFea,spcaTraFea,ldaTraFea,iprTraFea};
res.tesFea = {tesFea,pcaTesFea,spcaTesFea,ldaTesFea,iprTesFea};
end

function ScatterFeatures(traFea,tesFea,traCat,tesCat,method)
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
        ind = strcmp(uniCat(iCat),tesCat);
        scatter3(tesFea(ind,1),tesFea(ind,2),tesFea(ind,3),100,col{iCat},'Marker',mar{iCat},'MarkerFaceColor',col{iCat});
    end
end
title(['features transform: ' method]);
xlabel('Sepal length')
ylabel('Sepal width')

end

function cMat = TransformAndClassify(traFea,traCat,tesFea,tesCat,par)
cMat = confusionmat(nominal(tesCat),nominal(knnclassify(tesFea,traFea,traCat,par.nKNNs)));
end

function cMat = PCATransformAndClassify(traFea,traCat,tesFea,tesCat,par)
[traFeaNew,tesFeaNew] = PCAFeaturesTransform(traFea,tesFea);
cMat = confusionmat(nominal(tesCat),nominal(knnclassify(tesFeaNew,traFeaNew,traCat,par.nKNNs)));
end

function cMat = SPCATransformAndClassify(traFea,traCat,tesFea,tesCat,par)
[traFeaNew,tesFeaNew] = SPCAFeaturesTransform(traFea,tesFea,traCat);
cMat = confusionmat(nominal(tesCat),nominal(knnclassify(tesFeaNew,traFeaNew,traCat,par.nKNNs)));
end

function cMat = LDATransformAndClassify(traFea,traCat,tesFea,tesCat,par)
[traFeaNew,tesFeaNew] = LDAFeaturesTransform(traFea,tesFea,traCat);
cMat = confusionmat(nominal(tesCat),nominal(knnclassify(tesFeaNew,traFeaNew,traCat,par.nKNNs)));
end

function cMat = IPRTransformAndClassify(traFea,traCat,tesFea,tesCat,par)
[~,traFeaNew,tesFeaNew] = IPRLearnDisSub(traFea,tesFea,cellstr(traCat),par);
cMat = confusionmat(nominal(tesCat),nominal(knnclassify(tesFeaNew,traFeaNew,traCat,par.nKNNs)));
end