function plotorderparams(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

runTime=load('SCALARS/runTime.txt');

for i=1:numBlocks
    mindata=importdata(sprintf('SCALARS/minPhi%d.txt',i-1));
    meandata=importdata(sprintf('SCALARS/meanPhi%d.txt',i-1));
    maxdata=importdata(sprintf('SCALARS/maxPhi%d.txt',i-1));

    minPhi(i,:)=mindata';
    meanPhi(i,:)=meandata';
    maxPhi(i,:)=maxdata';
end

figure(3); 
minindex = length(runTime);
minindex = min(minindex, length(minPhi(1,:)));
t=runTime(1:minindex);
step=1;
tsub=runTime(1:step:minindex);

for i=1:numBlocks
    minsub(i,:)=minPhi(i,1:step:minindex);
    meansub(i,:)=meanPhi(i,1:step:minindex);
    maxsub(i,:)=maxPhi(i,1:step:minindex);

    op=(maxsub-minsub);
    rgb = hsv2rgb([(2/3)*(i-1)/(numBlocks-1) 1 1]);

    semilogy(tsub,op(i,:),'-','LineWidth',2,'Color',rgb); hold on;
end

hold off

xlabel('TIME')%,'FontWeight','bold');
ylabel('\phi_{max}-\phi_{min}')%,'FontWeight','bold');