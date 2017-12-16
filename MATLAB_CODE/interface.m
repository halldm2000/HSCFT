function interface(filename1,filename2,timeStep,rep)

temp = mgetfieldmpi3d(filename1,timeStep);
temp2d(:,:)= abs(temp(:,:,1));
field1 = repmat(temp2d,rep,rep);

temp = mgetfieldmpi3d(filename2,timeStep);
temp2d(:,:)= abs(temp(:,:,1));
field2 = repmat(temp2d,rep,rep);

runTime=load('SCALARS/runTime.txt');
[nx,ny]=size(field1);
%obtain append and write intervals
[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1})
aIndex = timeStep/appendInterval+1

h=surf(field1.*field2);
view(2);
axis tight;
pbaspect([ny nx 0.15*nx]);

material([.9 .1 0]);
colormap gray;
%caxis([0 1]);

lighting phong;
camlight left;
%shading flat;
shading interp;

camproj ortho;
