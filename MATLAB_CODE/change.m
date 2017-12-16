function cahnge(filename,time0,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= abs(temp(:,:,1));
field = repmat(temp2d,rep,rep);

temp = mgetfieldmpi3d(filename,time0);
temp2d(:,:)= abs(temp(:,:,1));
field0 = repmat(temp2d,rep,rep);

runTime=load('SCALARS/runTime.txt');
[nx,ny]=size(field);
%obtain append and write intervals
[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1})
aIndex = timeStep/appendInterval+1

h=surf(field-field0);
view(2);
axis tight;
pbaspect([ny nx 0.15*nx]);

material([1 0 0]);
colormap gray;
%caxis([0.1 .9]);
shading interp;
camproj ortho;
