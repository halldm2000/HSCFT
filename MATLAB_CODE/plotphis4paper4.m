function plotphis4paper4(timeStep)
figure(3); 

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})
runTime=load('SCALARS/runTime.txt');

symbol = {'ko','ks','kd','k^','k<','k>','kv'};
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
    step1=10;
    tsub1=runTime(1:step1:minindex);
    maxsub1=maxPhi(1:step1:minindex);
    minsub1=minPhi(1:step1:minindex);

    f=mean(meanPhi);
    step2=50;
    tsub2=runTime(1:step2:minindex);
    maxsub2=maxPhi(1:step2:minindex);
    minsub2=minPhi(1:step2:minindex);
    maxdiff = max(maxPhi(1:minindex)-minPhi(1:minindex));
    diff1 = maxsub1-minsub1;
    diff2 = maxsub2-minsub2;
    rgb = hsv2rgb([ (2/3)*(i-1)/(numBlocks-1), 1, 1]);
    %plot(tsub1,diff1,'-','Color','black','LineWidth',1.1); hold on
    plot(tsub2,diff2,char(symbol(i)),'LineWidth',1,'MarkerFaceColor','white','MarkerEdgeColor','black'); hold on
end

hold off
pbaspect([2 1 1])
xlabel('time','FontWeight','bold','FontSize',12);
ylabel('max(\phi)-min(\phi)','FontWeight','bold','FontSize',12);