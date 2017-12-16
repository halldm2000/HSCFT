function plotenergy(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

load SCALARS/runTime.txt
load SCALARS/S.txt
%load SCALARS/St0.txt;

for i=1:numBlocks
    data=importdata(sprintf('SCALARS/Sc%d.txt',i-1));
    Sc(i,:)=data';
end

figure(3); 
minindex = min( [ length(runTime) length(S) ]);

step=1; tsub=runTime(1:step:minindex);
plot(tsub, S(1:step:minindex),'-','LineWidth',1,'Color','k'); hold on;
%plot(tsub,St0(1:step:minindex),'-','LineWidth',2,'Color','k'); hold on;

for i=1:numBlocks
    rgb = hsv2rgb([(i-1)/(numBlocks-1) 1 1]);
    plot(tsub,Sc(i,:),'-','LineWidth',2,'Color',rgb); hold on;
end

hold off

xlabel('Time')%,'FontWeight','bold');
ylabel('Entropy Components')%,'FontWeight','bold');