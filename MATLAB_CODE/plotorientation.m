function [t,o]=plotorientation(field,timeSteps)
figure(4);

[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
t=timeSteps*writeInterval;
clear o;
for i=1:length(timeSteps)
    o(i)=orientation(field,timeSteps(i));
    figure(5); plot(t(1:i),o(1:i));
end

save orientation t o