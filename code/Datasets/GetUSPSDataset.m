function [fea,cat] = GetUSPSDataset(digits)

load usps_resampled

uniCat = strsplit(num2str(0:9),' ');
nCat = length(uniCat);
traTem = mod(find(train_labels==1),10);
tmpFea = [];
tmpCat = [];
for iCat=1:nCat
    traIdx = traTem==(str2double(uniCat{iCat}));
    tmpCat = [tmpCat; traTem(traIdx)];
    tmpFea = [tmpFea, train_patterns(:,traIdx)];
end
tmpCat(tmpCat==0) = 10;
tmpCat = tmpCat-1;
tmpFea = tmpFea';

tmpCat = cellstr(num2str(tmpCat));
cat = [];
fea = [];
if exist('digits','var')
    for iDig = 1:length(digits)
        idx = strcmp(num2str(digits(iDig)),tmpCat);
        cat = [cat; tmpCat(idx)];   
        fea = [fea; tmpFea(idx,:)];
    end
else
    fea = tmpFea;
    cat = tmpCat;
end

I = @(n) imagesc(reshape(fea(n,:),sqrt(size(fea,2)),sqrt(size(fea,2))));
clear train_labels train_patterns test_labels test_patterns traTem tesTem