function TestClassSubSpaRan

reset(RandStream.getGlobalStream) %reset random stream to reproduce results exactly

% add path and compile omp if not already done
addpath(genpath(pwd));
A = dir('SMALL/OMP_Private');
if length(strfind([A(:).name],'ompmex'))<3 %if the number of files containing 'ompmex' is less than 2 compile
    cd 'SMALL/OMP_Private'
    make
    cd ..
end

% parse dataset functions
datfun = {@GetToyExampleDataset,@GetFisherIrisDataset,@GetBalanceDataset,...
    @GetParkinsonsDataset,@GetSonarDataset};
datasets    = cell(length(datfun),1);
for i=1:length(datfun)
    s= functions(datfun{i});
    datasets{i} = regexprep(s.function,'Get(\w+)Dataset','$1');
end

% parameters
par.visu = false;

for iDat=1:length(datasets)
    fprintf('\nobtaining dataset %s... ',datasets{iDat});
    [fea,cat] = datfun{iDat}();                 %get features and categories
    if ~isfloat(fea), fea = double(fea); end    %transform features to float if needed
    fprintf('done\n');
    feaDim = size(fea,2);
    subSpaRanks = TrimEnd(ceil(linspace(1,feaDim,5)));
    for iSSR = 1:length(subSpaRanks)
        fprintf('testing transforms using subspaces of rank %d... ',subSpaRanks(iSSR));
        par.subSpaRan = subSpaRanks(iSSR);
        res(iDat,iSSR) = TestTransforms(fea,cat,par);
        fprintf('done\n');
    end
    mcrs = [res(iDat,:).mcr]';
    figure, PlotMCRs(mcrs(:,2:end),subSpaRanks,datasets{iDat});
end

function PlotMCRs(mcrs,subSpaRanks,dataset)
methodNames = {'PCA','SPCA','LDA','IPR'};
plot(subSpaRanks,mcrs,'o-','MarkerSize',10)
legend(methodNames);
xlabel('Sub-space rank');
ylabel('misclassification ratio');
title(sprintf('Classification of %s dataset',dataset));
grid on

function v = TrimEnd(v)
% trim the end of a linspace vector if it contains equal values
ind = find(v==max(v));
if length(ind) > 1
    v(ind(2:end)) = [];
end

