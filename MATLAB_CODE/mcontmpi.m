function mcontmpi(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);
runTime=load('SCALARS/runTime.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})
numLiquids = numBlocks-numSolids;

%obtain append and write intervals
[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));

writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});
aIndex = timeStep/appendInterval +1

cellval  = rvals( find(strcmp('L',rlabels)) );
L=str2num(cellval{1});

[ny,nx]=size(field)
lengthPerGP = L/(nx-1);
[x,y]=meshgrid(0:nx-1,0:ny-1);
[c,h]=contourf(x*lengthPerGP,y*lengthPerGP,1-field,1); hold on

pbaspect([nx ny .25]);
colormap gray

for i=numLiquids+1:1:numBlocks
    temp = mgetfieldmpi3d( sprintf('phi%d',i-1) ,timeStep);
    temp2d(:,:)= temp(:,:,1);
    wall = repmat(temp2d,rep);
    %[c,h]=contour(wall,[.7],'k');
    %set(h,'LineWidth',1.5);
end
hold off