function multisurf(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

timeScale =1;
vScale = 1;
marker={'.','o','*','s'}

rmsV=load('SCALARS/rmsV.txt');
rmsVTube=load('SCALARS/rmsVTube.txt');
runTime=load('SCALARS/runTime.txt');

for i=1:numBlocks
    data=importdata(sprintf('SCALARS/rmsW%d.txt',i-1));
    rmsW(i,:)=data';
end

figure(2); 
minindex = length(runTime);
minindex = min(minindex, length(rmsW(1,:)));
minindex = min(minindex, length(rmsV));

step=100;
tsub    = timeScale*runTime(1:step:minindex);
vsub    = vScale*rmsV(1:step:minindex);
vtsub   = vScale*rmsVTube(1:step:minindex);

[tsort,isort] = sort(tsub);
[t,tindex] = sort(runTime);

plot(t,rmsV(tindex),'k-','LineWidth',1.5); hold on
plot(t,rmsVTube(tindex),'k--','LineWidth',1.5);

for i=1:numBlocks
    wsub(i,:)=vScale*rmsW(i,1:step:minindex);
    if (numBlocks>1) rgb = hsv2rgb([(2/3)*(i-1)/(numBlocks-1) 1 1]);
    else rgb=[1 0 0];
    end;
    plot(t,rmsW(i,tindex),'-','LineWidth',1,'Color',rgb);
    plot(tsort,wsub(i,isort),char(marker(i)));
end

maxV = max(rmsV)
maxVIndex = find(rmsV==maxV);
maxVTime = runTime(maxVIndex);
pbaspect([2 1 1]);
xlabel('time');
