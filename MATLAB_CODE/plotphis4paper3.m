function plotphis(timeStep)
figure(3); 

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})
runTime=load('SCALARS/runTime.txt');

symbol = {'k-o','k-s','k-d','k-^','k-<','k->','k-v'};
labels = {'A','B','C','D','E'};

for i=1:numBlocks
    mindata=importdata(sprintf('SCALARS/minPhi%d.txt',i-1));
    meandata=importdata(sprintf('SCALARS/meanPhi%d.txt',i-1));
    maxdata=importdata(sprintf('SCALARS/maxPhi%d.txt',i-1));
    minPhi=mindata';
    meanPhi=meandata';
    maxPhi=maxdata';
    
    minindex = length(runTime);
    minindex = min(minindex, length(maxPhi));
    t=runTime(1:minindex);
    
    step1=50;
    tsub1=runTime(1:step1:minindex);
    maxsub1=maxPhi(1:step1:minindex);
    meansub1=meanPhi(1:step1:minindex);
    
    %semilogy(tsub1,(maxsub1-meansub1),char(symbol(i)),'LineWidth',1.3); hold on
    plot(tsub1,(maxsub1-meansub1),char(symbol(i)),'LineWidth',1.3); hold on
end

hold off
pbaspect([1.5 1 1])
xlabel('Time','FontWeight','bold','FontSize',12);
ylabel('Max(\phi)-Mean(\phi)','FontWeight','bold','FontSize',12);
legend(labels(1:numBlocks),4); legend boxoff;