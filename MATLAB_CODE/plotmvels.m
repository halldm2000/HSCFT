function multisurf(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

load SCALARS/rmsV.txt;
load SCALARS/rmsVTube.txt;
load SCALARS/runTime.txt;

for i=1:numBlocks
    data=importdata(sprintf('SCALARS/rmsW%d.txt',i-1));
    rmsW(i,:)=data';
end

figure(2); 
minindex = length(runTime);
minindex = min(minindex, length(rmsW(1,:)));
minindex = min(minindex, length(rmsV));

t=runTime(1:minindex);
step=1;
tsub=runTime(1:step:minindex);
vsub=rmsV(1:step:minindex);
vTubeSub=rmsVTube(1:step:minindex);
semilogy(tsub,vsub,'k--','LineWidth',2);
hold on

semilogy(tsub,vTubeSub,'k-','LineWidth',2);

for i=1:numBlocks
    wsub(i,:)=rmsW(i,1:step:minindex);
    rgb = hsv2rgb([(2/3)*(i-1)/(numBlocks-1) 1 1]);
    semilogy(tsub,wsub(i,:),'-','LineWidth',2,'Color',rgb);
end

hold off

xlabel('TIME')%,'FontWeight','bold');
ylabel('VELOCITIES')%,'FontWeight','bold');

maxV = max(rmsV)
maxVIndex = find(rmsV==maxV);
maxVTime = runTime(maxVIndex);