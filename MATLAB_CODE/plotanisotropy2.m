function [a,t,o]=plotanisotropy2(field,timeSteps)

[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
t=timeSteps*writeInterval;
clear a; clear o;
for i=1:length(timeSteps)
    [a(i),o(i)]=anisotropy2(field,timeSteps(i));
    figure(4); plot(t(1:i),a(1:i));
    figure(5); plot(t(1:i),o(1:i));
end

save orientation a t o