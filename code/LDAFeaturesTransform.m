function [newTraFea, newTesFea] = LDAFeaturesTransform(traFea,tesFea,traCat)
if ~nargin, unitTest; return; end

ldaObj = ClassificationDiscriminant.fit(traFea,traCat);           %fit a LDA object
Sb = ldaObj.BetweenSigma;                   %Between classes scatter matrix
Sw = ldaObj.Sigma;                          %Within classes scatter matrix
St = Sb + Sw;                               %Total scatter matrix

J = @(W) trace(pinv(W'*St*W)*(W'*Sb*W));    %as a function of projection matrix

[D,X] = eig(pinv(St)*Sb);                   %compute evd
idx = diag(real(X))>1e-8;                   %find non-zero eigenvalues
W = real(D(:, idx));                        %select corresponding eigenvectors
W = real(D(:,1:length(unique(traCat))-1));
% assert(size(W,2)==(length(unique(traCat))-1));

newTraFea = (W*pinv(W)*traFea')';           %project data onto subspace spanned by eigenvectors
newTesFea = (W*pinv(W)*tesFea')';


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

[newTraFea,newTesFea] = LDAFeaturesTransform(traFea,tesFea,traCat);

figure
gscatter(newTraFea(:,1), newTraFea(:,2), traCat,'rgb','osd');
xlabel('Sepal length');
ylabel('Sepal width');

J1 = @(v) trace(pinv(ClassificationDiscriminant.fit(v,traCat).Sigma)*ClassificationDiscriminant.fit(v,traCat).BetweenSigma);
J2 = @(v) trace(ClassificationDiscriminant.fit(v,traCat).BetweenSigma)/trace(ClassificationDiscriminant.fit(v,traCat).Sigma);
J3 = @(v) log(abs(det(ClassificationDiscriminant.fit(v,traCat).BetweenSigma))) - ClassificationDiscriminant.fit(v,traCat).LogDetSigma;

figure, subplot(3,1,1)
plot([J1(traFea) J1(newTraFea)],'*k')
xlim([.5 2.5])

subplot(3,1,2)
plot([J2(traFea) J2(newTraFea)],'*k'), hold on
xlim([.5 2.5])

subplot(3,1,3)
plot([J3(traFea) J3(newTraFea)],'*k'), hold on
xlim([.5 2.5])