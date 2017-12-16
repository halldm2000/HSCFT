function plotevfs4paper(plotnum);

figure(3);
set(0,'DefaultAxesColorOrder',[0 0 0]);
style={'.','o','*','s'}
    
[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

load SCALARS/runTime.txt;
data=importdata(sprintf('SCALARS/evf0.txt',i-1))';
L = min(length(runTime),length(data));


step=200; t1=runTime(1:step:L); evf1=data(1:step:L);
step=100; t2=runTime(1:step:L); evf2=data(1:step:L);

h=plot(t2,evf2,'k-'); hold on
h=plot(t1,evf1,char(style(plotnum)),'MarkerFaceColor','w'); hold on

xlabel('Time')
ylabel('Excess Volume Fraction')
