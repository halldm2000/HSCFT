function plotlengthscale(plotnum);

figure(3);
style={'ko','k.','ks','k*'}
  
[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

load SCALARS/runTime.txt;
load SCALARS/Rmean0.txt;
L = min(length(runTime),length(Rmean0))

step=10; t1=runTime(1:step:L); l1=Rmean0(1:step:L);
step=1000; t2=runTime(1:step:L); l2=Rmean0(1:step:L);

h=plot(t1,l1,'k-'); hold on
h=plot(t2,l2,char(style(plotnum)),'MarkerFaceColor','w'); hold on
%h=plot(t2,l2,char(style(plotnum)) ); hold on

xlabel('Time')
ylabel('Characteristic Length Scale')
pbaspect([1 1 1]);
axis tight

maxLength=max(max(Rmean0))