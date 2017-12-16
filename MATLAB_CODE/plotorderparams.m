function plotorderparams(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

runTime=load('SCALARS/runTime.txt');
figure(3); 
minindex = length(runTime);
for i=1:numBlocks
    mindata=importdata(sprintf('SCALARS/minPhi%d.txt',i-1));
    meandata=importdata(sprintf('SCALARS/meanPhi%d.txt',i-1));
    maxdata=importdata(sprintf('SCALARS/maxPhi%d.txt',i-1));

    minPhi(i,:)=mindata';
    meanPhi(i,:)=meandata';
    maxPhi(i,:)=maxdata';
    minindex = min(minindex, length(minPhi(i,:)));
end

t=runTime(1:minindex);
step=1;
tsub=runTime(1:step:minindex);

for i=1:numBlocks
    minsub(i,:)=minPhi(i,1:step:minindex);
    meansub(i,:)=meanPhi(i,1:step:minindex);
    maxsub(i,:)=maxPhi(i,1:step:minindex);
    op(i,:)=(maxsub(i,:)-minsub(i,:));
end

step1=1;
tsub1=runTime(1:step1:minindex);

step2=round(minindex/10);
tsub2=runTime(1:step2:minindex);

msize = 10;
mwidth =1.5;
lwidth =1.5;
plot(tsub2,op(1,1:step2:minindex),'ok','LineWidth',mwidth,'MarkerFaceColor','r','MarkerSize',msize); hold on;
if (numBlocks>1 )plot(tsub2,op(2,1:step2:minindex),'sk','LineWidth',mwidth,'MarkerFaceColor','g','MarkerSize',msize); hold on; end;
if (numBlocks>2) plot(tsub2,op(3,1:step2:minindex),'^k','LineWidth',mwidth,'MarkerFaceColor','b','MarkerSize',msize); hold on; end;
if (numBlocks>3) plot(tsub2,op(4,1:step2:minindex),'xk','LineWidth',mwidth,'MarkerFaceColor','w','MarkerSize',msize); hold on; end;
if (numBlocks>4) plot(tsub2,op(5,1:step2:minindex),'+k','LineWidth',mwidth,'MarkerFaceColor','k','MarkerSize',msize); hold on; end;
if (numBlocks>5) plot(tsub2,op(6,1:step2:minindex),'.k','LineWidth',mwidth,'MarkerFaceColor','k','MarkerSize',msize); hold on; end;

plot(tsub1,op(:,1:step1:minindex),'-','LineWidth',lwidth,'Color','k'); hold on;
pbaspect([3 1 1]);
hold off
axis tight;
%ylim([0 1]);
if (numBlocks==2) legend('A','B'); end;
if (numBlocks==3) legend('A','B','C'); end;
if (numBlocks==4) legend('A','B','C','D'); end;
if (numBlocks==5) legend('A','B','C','D','E'); end;
if (numBlocks==6) legend('A','B','C','D','E','F'); end;



xlabel('t/\tau_d')%,'FontWeight','bold');
ylabel('\phi_{max}-\phi_{min}')%,'FontWeight','bold');