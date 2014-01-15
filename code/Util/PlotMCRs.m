function PlotMCRs(mcrs,subSpaRanks,nUniCat,dataset)
methodNames = {'none','LDA','PCA','S-PCA','S-IPR'};
plot(subSpaRanks(end),mean(mcrs(:,1)),'ks','MarkerSize',10)
hold on
plot(nUniCat-1,mean(mcrs(:,4)),'rs','MarkerSize',10)
plot(subSpaRanks,mcrs(:,[1,2,5]),'--o','MarkerSize',10)
legend(methodNames);
xlabel('Sub-space rank');
ylabel('misclassification ratio');
title(sprintf('Classification of %s dataset',dataset));
randomChoice = 1/nUniCat;
ylim([0,max([randomChoice, max(mcrs(:))])*1.2]);
plot(get(gca,'XLim'),[randomChoice,randomChoice],'k-');
text(min(get(gca,'xlim')),randomChoice*1.1,'chance level','Margin',10);
set(gca,'YGrid','on','XGrid','off',...
    'XTick',rmd([subSpaRanks(1:end-1), subSpaRanks(end)-1, subSpaRanks(end)]));
a = fileparts(which(mfilename));
savegraph([a filesep dataset])
saveas(gcf,[dataset '.fig'])