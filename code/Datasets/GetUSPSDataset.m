function [fea,cat] = GetUSPSDataset

load usps_resampled

uniCat = strsplit(num2str(0:9),' ');
nCat = length(uniCat);
traTem = mod(find(train_labels==1),10);
fea = [];
cat = [];
for iCat=1:nCat
    traIdx = traTem==(str2double(uniCat{iCat}));
    cat = [cat; traTem(traIdx)];
    fea = [fea, train_patterns(:,traIdx)];
end
cat(cat==0) = 10;
cat = cat-1;
fea = fea';
I = @(n) imagesc(reshape(fea(n,:),sqrt(size(fea,2)),sqrt(size(fea,2))));
clear train_labels train_patterns test_labels test_patterns traTem tesTem