function msurfmpi(filename,timeStep,rep)

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
h=surf(field);

view(3);
axis tight;
pbaspect([rep(2)*ny rep(1)*nx .25*nx]);
material([0.4 .8 .3]);
lighting phong;
set(h,'FaceColor',[.9 .9 .9]);
%shading interp; colormap jet; caxis([0 1]);
camlight left;
camproj perspective;
zlim([0 1]);
%view([90 0]);