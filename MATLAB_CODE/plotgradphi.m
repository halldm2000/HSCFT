function plotgradphi(timeStep,rep)

[labels,vals]=textread('DomainSize_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

load runTime.txt

for i=1:numBlocks
    data=importdata(sprintf('rmsGradPhi%d.txt',i-1));
    gradPhi(i,:)=data';
end

figure(3); 
minindex = length(runTime);
minindex = min(minindex, length(gradPhi(1,:)));

t=runTime(1:minindex);
step=1;
tsub=runTime(1:step:minindex);

for i=1:numBlocks
    gsub(i,:)=gradPhi(i,1:step:minindex);
    rgb = hsv2rgb([(i-1)/numBlocks 1 1]);
    plot(tsub,gsub(i,:),'LineWidth',2,'Color',rgb);
    if (i==1) hold on
    end
end

hold off
for i=1:numBlocks
    maxGrad(i)=max(max(gradPhi(i,:)));
end
maxGrad
xlabel('TIME')
ylabel('PHI GRADIENTS')