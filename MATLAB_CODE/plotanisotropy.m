function [a,t]=plotanisotropy(field,timeSteps)

[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
t=timeSteps*writeInterval;
clear a;
for i=1:length(timeSteps)
    a(i)=anisotropy(field,timeSteps(i));
    figure(4); plot(t(1:i),a(1:i));
end