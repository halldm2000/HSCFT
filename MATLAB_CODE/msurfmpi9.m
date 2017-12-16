function msurfmpi(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);
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
h=surf(field);
%axis off
view(3);
axis tight;
pbaspect([rep(2)*ny rep(1)*nx .25*nx]);
lighting phong;
set(h,'FaceColor',[.99 .99 .99]);
light('Color',[.4 .4 .4],'Position',[-1 -1 1]);
camproj perspective;
light('Color',[0 0 .5],'Position',[0 0 1]);
light('Color',[0 .3  0],'Position',[-1 -1 0]);
%set(gcf,'Color','white')
%set(gca,'Color','white')

material([0 .99 .3]);

%zlim([0 1]);
%view([90 0]);