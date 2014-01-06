function [traFea, tesFea] = NormalizeFeatures(traFea,tesFea)
m = mean(traFea,1);
s = std(traFea,[],1);
for iFea=1:size(traFea,2)
    traFea(:,iFea) = (traFea(:,iFea)-m(iFea))./max(s(iFea),realmin);
end
if exist('tesFea','var')
    for iFea=1:size(tesFea,2)
        tesFea(:,iFea) = (tesFea(:,iFea)-m(iFea))./s(iFea);
    end
end