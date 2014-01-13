function hh = gscatter3(x,y,z,cat)
% 3d group scatter
uniCat = unique(cat);
nUniCat = length(uniCat);
for iUniCat = 1:nUniCat
    ind = find(strcmp(cat,uniCat(iUniCat)));
    hh = scatter3(x(ind,:),y(ind,:),z(ind,:));
    hold on
end