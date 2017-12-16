function plotScatteringFunction
load SCALARS/scatteringFcn.txt;
[nSteps,kBins]=size(scatteringFcn);

figure(4)
 
plot(scatteringFcn')
[sx,sy]=size(scatteringFcn')
x=0:kBins-1;
xm = ones(sy,1)*x
meanx = x*scatteringFcn'./sum(scatteringFcn')
%semilogy( (1.0./meanx)*xm,scatteringFcn')

loglog( (1.0./meanx)*xm,scatteringFcn')
axis tight;