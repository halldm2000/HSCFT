function plotphis2(timeStep)
figure(3); 

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})
runTime=load('SCALARS/runTime.txt');

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
    step=1;
    tsub=runTime(1:step:minindex);
    maxsub=maxPhi(1:step:minindex);
    meansub=meanPhi(1:step:minindex);
    diff = maxsub-meansub;
    maxdiff= max(diff);
    rgb = hsv2rgb([ (2/3)*(i-1)/(numBlocks-1), 1, 1]);
    %plot(tsub,diff/maxdiff,'-','LineWidth',1.2,'Color',rgb); hold on
     semilogy(tsub,diff-diff(1),'-','LineWidth',1.2,'Color',rgb); hold on
end

hold off

xlabel('time')%,'FontWeight','bold');
ylabel('Max(\phi_i)')%,'FontWeight','bold');