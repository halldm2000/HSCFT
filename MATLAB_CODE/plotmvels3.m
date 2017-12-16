function plotmvels3(timeStep,rep)

[labels,vals]=textread('DomainSize_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

[labels,vals]= textread('DIMENSIONAL_PARAMETERS.txt','%s\t%s');
cellval     = vals( find(strcmp('CONVECTION_TIME',labels)) );
timeScale   = str2num(cellval{1})
cellval     = vals( find(strcmp('CONVECTION_SPEED',labels)) );
vScale      = str2num(cellval{1})*10^7

load rmsV.txt
load rmsVTube.txt
load runTime.txt

for i=1:numBlocks
    data=importdata(sprintf('rmsW%d.txt',i-1));
    rmsW(i,:)=data';
end

figure(2); 
minindex = length(runTime);
minindex = min(minindex, length(rmsW(1,:)));
minindex = min(minindex, length(rmsV));

step=5;
tsub    = timeScale*runTime(1:step:minindex);
vsub    = vScale*rmsV(1:step:minindex);
vtsub   = vScale*rmsVTube(1:step:minindex);

plot(tsub,vsub,'k-o','MarkerSize',3,'LineWidth',1);
hold on

for i=1:numBlocks
    wsub(i,:)=vScale*rmsW(i,1:step:minindex);
    rgb     = hsv2rgb([(i-1)/(numBlocks) 1 1]);
    plot(tsub,wsub(i,:),'-s','MarkerSize',3,'Color',rgb);
end
axis tight;
hold off;

xlabel('sec')%,'FontWeight','bold');
ylabel('nm/sec')%,'FontWeight','bold');

maxV = max(rmsV)
maxVIndex = find(rmsV==maxV);
maxVTime = runTime(maxVIndex);

for i=1:numBlocks
    maxW(i)=max(rmsW(i,:));
    maxWIndex=find(rmsW(i,:)==maxW(i));
    maxWTime(i)=runTime(maxWIndex);
end
maxWTime