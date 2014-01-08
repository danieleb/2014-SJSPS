function TestClassSubSpaRan
subSpaRanks = 1:4;
datasets    = {'coloncancer'};
datfun = @GetColonCancerDataset;
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
