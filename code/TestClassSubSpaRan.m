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
datfun = {@GetFisherIrisDataset,@GetBalanceDataset,...
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
    subSpaRanks = rmd(ceil(linspace(1,feaDim,22)));
    for iSSR = 1:length(subSpaRanks)
        fprintf('testing transforms using subspaces of rank %d... ',subSpaRanks(iSSR));
        par.subSpaRan = subSpaRanks(iSSR);
        res(iDat,iSSR) = TestTransforms(fea,cat,par);
        fprintf('done\n');
    end
    mcrs = [res(iDat,:).mcr]';
    figure, PlotMCRs(mcrs,subSpaRanks,length(unique(cat)),datasets{iDat});
end

