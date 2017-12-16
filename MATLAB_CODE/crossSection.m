function crossSection(filename,timeStep,index,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);
runTime=load('SCALARS/runTime.txt');

%[labels,vals]= textread('DIMENSIONAL_PARAMETERS.txt','%s\t%s');
%cellval      = vals( find(strcmp('LENGTH',labels)) );
%lengthScale    = str2num(cellval{1})*10^7

%obtain append and write intervals
[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});
aIndex = timeStep/appendInterval +1

[nx,ny]=size(field);
h=plot(field(:,index),'k','LineWidth',2);


ylim([0 1]);
