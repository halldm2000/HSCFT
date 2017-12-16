function isophis(field,t,rep)
filename='phi0'

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});
cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1});
numLiquids = numBlocks-numSolids;

load SCALARS/L.txt;

temp = mgetfieldmpi3d(field,t); 
temp = repmat(temp,rep);
phi(:,:,:) = temp;

[nr,nc,nd]=size(temp);
n=size(temp); ny=nr; nx=nc; nz=nd;

figure(1); cla;
meanphi = mean(mean(mean(phi)))

rgb=[.99 .99 .99];
rgb2=[.4 .4 .4];

% plot the isosurfaces for all fields included in the array "which"
vx=[0:nx-1]*rep(2)*L(t+1,1)/(nx-1);
vy=[0:ny-1]*rep(1)*L(t+1,2)/(ny-1);
vz=[0:nz-1]*rep(3)*L(t+1,3)/(nz-1);
[x,y,z]=meshgrid(vx,vy,vz);

p1 = patch(isosurface(x,y,z,phi,meanphi,phi),'FaceColor',rgb,'EdgeColor','none');
p2 = patch(isocaps(x,y,z,phi,meanphi),'FaceColor',rgb2,'EdgeColor','none');
isonormals(x,y,z,phi,p1);
xlabel('x'); ylabel('y'); zlabel('z');

box on;
axis equal tight;
xlim([0 rep(2)*L(t+1,1)]);
ylim([0 rep(1)*L(t+1,2)]);
zlim([0 rep(3)*L(t+1,3)]);

view(3);
lighting phong;
material([0.1 1 .4]);
camproj perspective;
camlight(-40 ,40);

timestep=t
L(t+1,:)