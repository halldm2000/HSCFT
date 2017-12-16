function plotevf(timeStep,rep)

[labels,vals]=textread('DomainSize_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

load runTime.txt

for i=1:numBlocks
    data=importdata(sprintf('evf%d.txt',i-1));
    gradPhi(i,:)=data';
end

figure(3); 
minindex = length(runTime);
minindex = min(minindex, length(gradPhi(1,:)));

t=runTime(1:minindex);
step1=30;
tsub1=runTime(1:step1:minindex);
gsub1(:,:)=gradPhi(:,1:step1:minindex);

step2=100;
tsub2=runTime(1:step2:minindex);
gsub2(:,:)=gradPhi(:,1:step2:minindex);

rgb = hsv2rgb([(2/3)*0/(numBlocks-1) 1 1]);
plot(tsub2,gsub2(1,:),'^k','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*1/(numBlocks-1) 1 1]);
plot(tsub2,gsub2(2,:),'ok','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*2/(numBlocks-1) 1 1]);
plot(tsub2,gsub2(3,:),'sk','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*3/(numBlocks-1) 1 1]);
plot(tsub2,gsub2(4,:),'dk','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*4/(numBlocks-1) 1 1]);
plot(tsub2,gsub2(5,:),'vk','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*0/(numBlocks-1) 1 1]);
plot(tsub1,gsub1(1,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*1/(numBlocks-1) 1 1]);
plot(tsub1,gsub1(2,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*2/(numBlocks-1) 1 1]);
plot(tsub1,gsub1(3,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*3/(numBlocks-1) 1 1]);
plot(tsub1,gsub1(4,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*4/(numBlocks-1) 1 1]);
plot(tsub1,gsub1(5,:),'-','LineWidth',2,'Color',rgb); hold on;

hold off
pbaspect([2 1 1]);
axis tight;
ylim([0.25 .65]);

xlabel('Time');
ylabel('Excess Area Fraction');
legend('block A','block B','block C','block D','block E');