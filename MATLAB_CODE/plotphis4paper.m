function plotphis(timeStep)
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
    step1=100;
    tsub1=runTime(1:step1:minindex);
    maxsub1=maxPhi(1:step1:minindex);
    
    step2=200;
    tsub2=runTime(1:step2:minindex);
    maxsub2=maxPhi(1:step2:minindex);
    rgb = hsv2rgb([ (2/3)*(i-1)/(numBlocks-1), 1, 1]);
    plot(tsub1,maxsub1,'k-','LineWidth',1.2); hold on
    plot(tsub2,maxsub2,char(symbol(i)),'LineWidth',1,'MarkerFaceColor',rgb,'MarkerSize',8); hold on
end

hold off
pbaspect([1.5 1 1])
xlabel('time','FontWeight','bold','FontSize',12);
ylabel('max(\phi_i)','FontWeight','bold','FontSize',12);