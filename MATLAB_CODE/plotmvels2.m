function Wi=plotmvels2(timeStep,rep)

figure(2); 

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})

numLiquids = numBlocks-numSolids;

timeScale =1;
vScale = 1;

rmsV=load('SCALARS/rmsV.txt');
rmsVTube=load('SCALARS/rmsVTube.txt');
runTime=load('SCALARS/runTime.txt');

minindex = length(runTime);
minindex = min(minindex, length(rmsV));
minindex = min(minindex, length(rmsVTube));

step=1;
tsub    = timeScale*runTime(1:step:minindex);
vsub    = vScale*rmsV(1:step:minindex);
vtsub   = vScale*rmsVTube(1:step:minindex);

[tsort,isort] = sort(tsub);
plot(tsort,vsub(isort),'k--','LineWidth',1.2); hold on
plot(tsort,vsub(isort),'k-','LineWidth',1.2); hold on;

for i=1:numLiquids
    data=importdata(sprintf('SCALARS/rmsW%d.txt',i-1));
    rmsW=data';
    minindex = min(minindex, length(rmsW));
    wsub=vScale*rmsW(1:step:minindex);
    if (numBlocks>1) rgb = hsv2rgb([(2/3)*(i-1)/((numLiquids-1)+(numLiquids==1)) 1 1]);
    else rgb=[1 0 0];
    end;
    plot(tsort,wsub(isort),'-','LineWidth',1.2,'Color',rgb);
end
hold off

maxV = max(rmsV)
tau=1;
L=load('SCALARS/L.txt');
Wi=maxV*tau/L(1,1)
maxVIndex = find(rmsV==maxV);
maxVTime = runTime(maxVIndex);
pbaspect([2 1 1]);
