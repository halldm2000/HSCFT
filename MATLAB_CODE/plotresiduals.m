function plotomegas

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

runTime=load('SCALARS/runTime.txt');

for i=1:numBlocks
    odata=importdata(sprintf('SCALARS/omegaErr%d.txt',i-1));
    omegaErr(i,:)=odata';
end
pdata=importdata(sprintf('SCALARS/rmsVDiv.txt',i-1));
vdata=importdata(sprintf('SCALARS/rmsVErr.txt',i-1));
ddata=importdata(sprintf('SCALARS/rmsdVdL.txt',i-1));
mvdata=importdata(sprintf('SCALARS/meanV.txt',i-1));

pErr=pdata';
vErr=vdata';
dErr = ddata';
mvErr=mvdata';

figure(3); 
minindex = length(runTime);
for i=1:numBlocks 
    minindex = min(minindex, length(omegaErr(i,:)));
end
minindex = min(minindex, length(pErr));

t=runTime(1:minindex);
step=1;
tsub=runTime(1:step:minindex);
pSub=pErr(1:step:minindex);
vSub=vErr(1:step:minindex);
dSub=dErr(1:step:minindex);
mvSub=mvErr(1:step:minindex);

for i=1:numBlocks
    omegaSub(i,:)=omegaErr(i,1:step:minindex);
    rgb = hsv2rgb([(i-1)/(numBlocks) 1 1]);
    semilogy(tsub,omegaSub(i,:),'-','LineWidth',1,'Color',rgb);
    if (i==1) hold on
    end
end
semilogy(tsub,pSub,'LineWidth',2,'Color','b');
semilogy(tsub,vSub,'LineWidth',2,'Color','g');
semilogy(tsub,dSub,'LineWidth',2,'Color','m');
semilogy(tsub,mvSub,'LineWidth',2,'Color','y');

hold off

xlabel('TIME')%,'FontWeight','bold');
ylabel('RESIDUAL ERRORS')%,'FontWeight','bold');
