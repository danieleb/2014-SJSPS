function [fea,cat] = GetColonCancerDataset()
% Gets the colon cancer dataset from the web
url = 'http://genomics-pubs.princeton.edu/oncology/affydata/I2000.html';
s = urlread(url);
s = regexprep(s,'<.*?>','');
format = repmat('%f',1,62);
C = textscan(s(26:end),format,'CollectOutput',true);
fea = C{1}';

catUrl = 'http://genomics-pubs.princeton.edu/oncology/affydata/tissues.html';
s = urlread(catUrl);
s = regexprep(s,'<.*?>','');
format = '%f';
C = textscan(s(26:end),format);
cat = cell(size(C{1},1),1);
for i=1:length(cat)
    if C{1}(i)<0
        cat{i} = 'tumor';
    else
        cat{i} = 'normal';
    end
end
