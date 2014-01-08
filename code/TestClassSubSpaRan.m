function TestClassSubSpaRan
% add path and compile omp if not already done
addpath(genpath(pwd));
A = dir('SMALL/OMP_Private');
if length(strfind([A(:).name],'ompmex'))<3 %if the number of files containing 'ompmex' is less than 2 compile
    cd 'SMALL/OMP_Private'
    make
    cd ..
end

subSpaRanks = 1:4;
datasets    = {'fisher'};
datfun = @GetFisherIrisDataset;
par.visu = false;
for iDat=1:length(datasets)
    fprintf('\nobtaining dataset %s... ',datasets{iDat});
    [fea,cat] = datfun();
    if ~isfloat(fea), fea = double(fea); end     %transform features to float if needed
    fprintf('done\n');
    for iSSR = 1:length(subSpaRanks)
        fprintf('testing transforms using subspaces of rank %d... ',subSpaRanks(iSSR));
        par.subSpaRan = iSSR;
        res(iSSR) = TestTransforms(fea,cat,par);
        fprintf('done\n');
    end
    mcrs = [res(:).mcr]';
    figure, PlotMCRs(mcrs(:,2:end));
end

function PlotMCRs(mcrs)
methodNames = {'PCA','SPCA','LDA','IPR'};
plot(mcrs,'o-','MarkerSize',10)
legend(methodNames);
xlabel('Sub-space rank');
ylabel('misclassification ratio');
grid on
