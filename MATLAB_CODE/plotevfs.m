function plotevf(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

runTime=load('SCALARS/runTime.txt');
figure(3)
for i=1:numBlocks
    data=importdata(sprintf('SCALARS/evf%d.txt',i-1))';
    rgb = hsv2rgb([(2/3)*(i-1)/(numBlocks-1) 1 1]); hold on
    L = min( length(data), length(runTime));
    plot(runTime(1:L),data(1:L),'-','LineWidth',1.5,'Color',rgb);
end

%hold off
ylim([0.0 1]);

xlabel('Time')%,'FontWeight','bold');
ylabel('Excess Volume Fraction')%,'FontWeight','bold');