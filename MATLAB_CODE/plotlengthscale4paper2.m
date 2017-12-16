function plotlengthscale(plotnum);

figure(3);
style={'.','o','*','s'}
  
[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})
for i=1:numBlocks
    data=importdata(sprintf('SCALARS/rmsW%d.txt',i-1));
    rmsW(i,:)=data';
end
load SCALARS/runTime.txt;
load SCALARS/Rmean0.txt;
peakV = max(max(rmsW))
peakTime = runTime(find(rmsW==peakV))

L = min(length(runTime),length(Rmean0))

step=10; t1=runTime(1:step:L); l1=Rmean0(1:step:L);
step=300; t2=runTime(1:step:L); l2=Rmean0(1:step:L);

h=loglog(t1,l1,'-'); hold on
h=semilogx(t2,l2,char(style(plotnum)),'MarkerFaceColor','w'); hold on

xlabel('Time')
ylabel('Characteristic Length Scale')
pbaspect([2 1 1]);
axis tight