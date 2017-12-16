function plotorderparams(timeStep,rep)

[labels,vals]=textread('DomainSize_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

load runTime.txt

for i=1:numBlocks
    mindata=importdata(sprintf('minPhi%d.txt',i-1));
    meandata=importdata(sprintf('meanPhi%d.txt',i-1));
    maxdata=importdata(sprintf('maxPhi%d.txt',i-1));

    minPhi(i,:)=mindata';
    meanPhi(i,:)=meandata';
    maxPhi(i,:)=maxdata';
end

figure(3); 
minindex = length(runTime);
minindex = min(minindex, length(minPhi(1,:)));
t=runTime(1:minindex);

bigStep=200;
smallStep=10;
tsub1=runTime(1:bigStep:minindex);
tsub2=runTime(1:smallStep:minindex);

op=maxPhi-minPhi;
opSub1(:,:)=op(:,1:bigStep:minindex);
opSub2(:,:)=op(:,1:smallStep:minindex);

rgb = hsv2rgb([(2/3)*0/(numBlocks-1) 1 1]);
plot(tsub1,opSub1(1,:),'^k','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*1/(numBlocks-1) 1 1]);
plot(tsub1,opSub1(2,:),'ok','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*2/(numBlocks-1) 1 1]);
plot(tsub1,opSub1(3,:),'sk','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*3/(numBlocks-1) 1 1]);
plot(tsub1,opSub1(4,:),'dk','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*4/(numBlocks-1) 1 1]);
plot(tsub1,opSub1(5,:),'vk','LineWidth',1.5,'MarkerFaceColor',rgb,'MarkerSize',7); hold on;

rgb = hsv2rgb([(2/3)*0/(numBlocks-1) 1 1]);
plot(tsub2,opSub2(1,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*1/(numBlocks-1) 1 1]);
plot(tsub2,opSub2(2,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*2/(numBlocks-1) 1 1]);
plot(tsub2,opSub2(3,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*3/(numBlocks-1) 1 1]);
plot(tsub2,opSub2(4,:),'-','LineWidth',2,'Color',rgb); hold on;

rgb = hsv2rgb([(2/3)*4/(numBlocks-1) 1 1]);
plot(tsub2,opSub2(5,:),'-','LineWidth',2,'Color',rgb); hold on;

hold off
axis tight;
pbaspect([2 1 1]);

xlabel('Time')
ylabel('\phi_{max}-\phi_{min}')
legend('block A','block B','block C','block D','block E');