function radius=frontradius(t,draw)
L = load('SCALARS/L.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1});
numLiquids = numBlocks-numSolids;

amplitude=-1;
amplitude2=-1;

for i=1:numLiquids
    % read the data
    data = mgetfieldmpi3d( sprintf('phi%d',i-1) ,t);
    phi(:,:)=data(:,:,1);
    
    % set surface amplitude to max concentration
    amplitude=max(amplitude,phi>0.501);
    amplitude2=max(amplitude2,phi);

end
rep=[1 1 1];
[nx,ny]=size(phi);
vx=[0:nx-1]*(L(t+1,2)/(nx-1));
vy=[0:ny-1]*(L(t+1,1)/(ny-1));
[x,y]=meshgrid(vy,vx);
d=sqrt( (x-mean(vx)).^2 + (y-mean(vy)).^2 ).*amplitude;
radius = max(max(d))

if(draw)
    figure(1); surf(d);
    shading interp; axis tight;
    pbaspect([ny nx 0.25*nx] )
    view(2);
    lighting phong;
    camlight;
    material dull;
end;