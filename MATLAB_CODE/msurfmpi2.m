function msurfmpi(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= abs(temp(:,:,1));
field = repmat(temp2d,rep,rep);
runTime=load('SCALARS/runTime.txt');
[nx,ny]=size(field);
%obtain append and write intervals
if exist('./STARTUP_FILES/RunParameters.txt') [rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
else [rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s'); end;writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1})
aIndex = timeStep/appendInterval+1

h=surf(field);
view(2);
axis tight;
%pbaspect([1 1 1]);
pbaspect([ny nx 0.15*nx]);

%material([.3 .7 .3]);
material dull;
%load redwhite; colormap(redwhite);
%load graywhite; colormap(graywhite);

caxis([-.01 .01]);
lighting phong;
%camlight left;
shading interp;

camproj ortho;
%titleString=sprintf('t=%.0f',runTime(aIndex));
%title(titleString);