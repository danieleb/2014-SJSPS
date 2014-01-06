function [cat, C] = ComputeAtomsCategories(X,cat)

if isnumeric(cat), cat = num2str(cat); end

nAto = size(X,1);
uniCat = unique(cat);
nCat = length(uniCat);
%compute a matrix whose (i,j)th element is the contribution of atom i to the sparse approximation of signals in the class j
C = zeros(nAto,nCat);
for iCat=1:nCat
    idx = strcmp(cellstr(cat),char(uniCat(iCat)));
    C(:,iCat) = sum(abs(X(:,idx)'))/sum(idx);
end

[~, cat] = max(C,[],2);
