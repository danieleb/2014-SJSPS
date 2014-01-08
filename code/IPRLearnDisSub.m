function [traSubSpa,newTraFea,newTesFea] = IPRLearnDisSub(traFea,tesFea,traCat,par)
% Learn discriminative subspaces using iterative projections and rotations
% (IPR) incoherent dictionaty learning algorithm.

if ~nargin, unitTest; return; end
%% Params and defaults
if ~exist('par','var') || isempty(par), par=struct; end
if isnumeric(traCat), traCat = num2str(traCat);     end

uniCat = unique(traCat);                %names of categories
nCat = length(uniCat);                  %number of categories
nDim = size(traFea,2);                  %dimension of features
toolbox   = 'SMALL_incoherentDL';		%dictionary learning toolbox

def.subSpaRan = 50;                     %rank of incoherent subspaces
def.nAto = nCat*nDim;                   %number of atoms
def.nIte = 1;                           %number of DL iterations

def.mu = sqrt((par.nAto-nDim)/(nDim*(par.nAto-1))); %mutual coherence
par = setdefaultoptions(par,def);       %set parameters

nAto  = par.nAto;                %number of atoms in the dictionary
perActAto = par.perActAto;   %percentage of active atoms
subSpaRan = par.subSpaRan;

nActiveAtoms = fix(nDim/100*perActAto); %number of active atoms


%% Create SMALL structures for dictionary learning
SMALL.Problem.b  = traFea';
SMALL.Problem.b1 = SMALL.Problem.b;
SMALL.Problem.cat = traCat;

% omp2 sparse representation solver
ompParam = struct('X',SMALL.Problem.b,'epsilon',1e-6,'maxatoms',nActiveAtoms); %parameters
solver	 = SMALL_init_solver('ompbox','omp2',ompParam,false); %solver structure

SMALL.DL.toolbox = toolbox;
SMALL.DL.name = 'ksvd';
SMALL.DL.profile = true;
SMALL.DL.param.data = SMALL.Problem.b;
SMALL.DL.param.Tdata = nActiveAtoms;
SMALL.DL.param.dictsize = nAto;
SMALL.DL.param.nIte = par.nIte;
SMALL.DL.param.memusage = 'high';
SMALL.DL.param.solver = solver;
SMALL.DL.param.decFcn = 'iprDis';
SMALL.DL.param.coherence = par.mu;
SMALL.DL.param.initdict = {};
SMALL.DL.param = setdefaultoptions(SMALL.DL.param,par);
SMALL.DL = SMALL_learn(SMALL.Problem,SMALL.DL);

D = SMALL.DL.D;                     %dictionary containing incoherent subspaces

%% Training phase
% Assign subspaces to categories
SMALL.Problem.A = D;
solver = SMALL_solve(SMALL.Problem,solver);
Xtra = solver.solution;

[q,C] = ComputeAtomsCategories(Xtra,traCat);  %atom's contributions to categories
for iCat=1:nCat     %for every category
    [~,idx] = sort(C(:,iCat),'descend');    %sort atoms according to contributions to category
    traSubSpa{iCat} = D(:,idx(1:subSpaRan));%select subset of atoms in the dictionary as subspace
    c=1;                                    %counter
    while rank(traSubSpa{iCat})<subSpaRan   %if and until rank of subspace is less than predefined rank
        traSubSpa{iCat} = [traSubSpa{iCat}, D(:,idx(subSpaRan+c))];  %add next atom
        c = c+1;
    end
end

% Project training data onto subspaces
Project = @(Phi,x) Phi*pinv(Phi)*x;         %projection operator
dis = nan(nCat);                            %matrix containing distance between subspaces and groups of observations
catIdx = cell(nCat,1);
newTraFea = nan(size(traFea));
for iCat=1:nCat
    for jCat=1:nCat
        catIdx{iCat} = strcmp(cellstr(traCat),char(uniCat(iCat)));     %indexes corresponding to iCat category
        tem = Project(traSubSpa{jCat},traFea(catIdx{iCat},:)')';       %Projection of signals of iCat catogory onto jCat subspace
        dis(iCat,jCat) = norm(traFea(catIdx{iCat},:)-tem,'fro');       %distance between data belonging to class iCat and their projection onto subspace jCat
    end
    [~, idx] = min(dis(iCat,:));
    newTraFea(catIdx{iCat},:) = Project(traSubSpa{idx},traFea(catIdx{iCat},:)')'; % new features are projection of old ones onto nearest subspace
end

%% Test phase
%Assign points to subspaces
newTesFea = nan(size(tesFea));
temp = nan(size(tesFea,2),nCat);
dis = nan(nCat,1);
for iTesFea=1:length(tesFea)        %for every data point in the test set
    for iCat=1:nCat                 %for every caegory
        temp(:,iCat) = Project(traSubSpa{iCat}, tesFea(iTesFea,:)');   %projection of datapoint onto subspace
        dis(iCat) = norm(tesFea(iTesFea,:)'-temp(:,iCat));    %distance between original datapoint and projection
    end
    [~, idx] = min(dis);                    %find the nearest subspace to the datapoint
    newTesFea(iTesFea,:) = Project(traSubSpa{idx},tesFea(iTesFea,:)')'; % new features are projection of old ones onto nearest subspace
end
    
end

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

newTraFea = IPRFeaturesTransform(traFea,tesFea,traCat);
size(newTraFea,2);
figure
if size(newTraFea,2)<2
    gscatter(newTraFea(:,1), newTraFea(:,1), traCat,'rgb','osd');
else
    gscatter(newTraFea(:,1), newTraFea(:,2), traCat,'rgb','osd');
end
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
end