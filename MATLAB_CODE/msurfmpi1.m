function msurfmpi(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);
runTime=load('SCALARS/runTime.txt');

%[labels,vals]= textread('DIMENSIONAL_PARAMETERS.txt','%s\t%s');
%cellval      = vals( find(strcmp('LENGTH',labels)) );
%lengthScale    = str2num(cellval{1})*10^7

%obtain append and write intervals
[rlabels,rvals]=textread('RUN_PARAMETERS/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});
aIndex = timeStep/appendInterval +1

[nx,ny]=size(field);
h=surf(field);

view(3);
axis tight;
pbaspect([rep(2) rep(1) .25]);
material([0.2 .8 .3]);
lighting phong;
shading interp
colormap autumn
axis off
%set(h,'FaceColor',[.8 .8 .8]);
camlight left;
camproj perspective;
%titleString=sprintf('step=%d t=%f',timeStep,runTime(aIndex));
%zlim([0 1.5]);