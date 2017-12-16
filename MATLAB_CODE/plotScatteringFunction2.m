function plotScatteringFunction(timesteps)

figure(4)

for i=1:length(timesteps)
    tstep = timesteps(i)
    filename = sprintf('scatteringFcn/scatteringFcn_%d.txt',timesteps(i));
    scatteringFcn(i,:) = load(filename);
    [nSteps,kBins]=size(scatteringFcn);
end
plot(scatteringFcn'); 
%x=1:kBins; loglog(x,scatteringFcn')

axis tight;
hold off;