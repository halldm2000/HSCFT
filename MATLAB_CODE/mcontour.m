function msurfmpi(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(1,:,:);
field = repmat(temp2d,rep,rep);
load runTime.txt;
load L.txt;

%obtain append and write intervals
[rlabels,rvals]=textread('RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1})
aIndex = timeStep/appendInterval+1

h=contourf(field,1);
axis tight;
pbaspect([L(aIndex,1) L(aIndex,2) 0.15*L(aIndex,1)]);
%camlight ;
titleString=sprintf('t=%.0f',runTime(aIndex));
title(titleString);