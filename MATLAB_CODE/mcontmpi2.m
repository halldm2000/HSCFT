function mcontmpi(filename,timeStep,rep,resolution,color)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);

%obtain append and write intervals
[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});
aIndex = timeStep/appendInterval +1

[labels,vals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
cellval  = vals( find(strcmp('L',labels)) );
L=str2num(cellval{1});
[ny,nx]=size(field)
lengthPerGP = L/(nx-1);
[x,y]=meshgrid(0:nx-1,0:ny-1);

[c,h]=contour(x*lengthPerGP,y*lengthPerGP,(1-field),1,'k'); hold on
set(h,'LineWidth',1,'Color',color);
axis equal;% tight
xlim([0 (nx-1)*lengthPerGP]);
ylim([0 (ny-1)*lengthPerGP]);
