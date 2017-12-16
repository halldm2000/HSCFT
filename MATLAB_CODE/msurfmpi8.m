function msurfmpi8(filename,timeStep,timeStep0,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);

temp = mgetfieldmpi3d(filename,timeStep0);
temp2d(:,:)= temp(:,:,1);
field0 = repmat(temp2d,rep);
runTime=load('SCALARS/runTime.txt');

%[labels,vals]= textread('DIMENSIONAL_PARAMETERS.txt','%s\t%s');
%cellval      = vals( find(strcmp('LENGTH',labels)) );
%lengthScale    = str2num(cellval{1})*10^7

%obtain append and write intervals
[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});
aIndex = timeStep/appendInterval +1

[nx,ny]=size(field);
h=surf(field-field0);

view(3);
axis tight;
pbaspect([rep(2)*ny rep(1)*nx .15*nx]);
material([0.4 .8 .3]);
lighting phong;
set(h,'FaceColor',[.8 .8 .8]);
camlight left;
camproj perspective;
%zlim([0 1]);