function AnalyzeResults(res)
% Get dimensions
nDatasets = size(res,1);

% for every dataset, get the combination that led to the best ipr result
ssrs = {'1','25%','50%','75%','100%'};
mus = {'incoherent','coherent'};
nns = {'1','5','many'};
for iData = 1:nDatasets
    r = squeeze(res(iData,:,:,:)); %select dataset
    [~, nNNs, nMus] = size(r);
    mcrs = [r(:).mcr];  %build matrix of mcrs
    iprMCRs = mcrs(end,:);
    [~,ind] = min(iprMCRs);
    [ssr,nnn,mu] = MapIndex(ind,nNNs,nMus);
    fprintf('\nBest ipr with %s-rank subspace, %s neighbors and %s dictionary',...
        ssrs{ssr},nns{nnn},mus{mu});
end

function [ssr,nnn,mu] = MapIndex(n,nNNs,nMus)
ssr = ceil(n/prod([nNNs,nMus])-1)+1;
n = n-prod([nNNs,nMus])*(ssr-1);
nnn = ceil(n/prod(nMus)-1)+1;
n = n - prod(nMus)*(nnn-1);
mu =n;
% fprintf('\n%d - %d %d %d',n,ssr,nnn,mu)