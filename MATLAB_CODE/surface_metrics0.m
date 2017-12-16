function surface_metrics(field,t,rep)

%load data
[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) ); numBlocks=str2num(cellval{1});
cellval  = vals( find(strcmp('numSolids',labels)) ); numSolids=str2num(cellval{1});
numLiquids = numBlocks-numSolids;
load SCALARS/L.txt;
temp = mgetfieldmpi3d(field,t); temp = repmat(temp,rep); phi(:,:,:) = temp;
[nr,nc,nd]=size(temp);n=size(temp); ny=nr; nx=nc; nz=nd;

% construct an isosurface
figure(1); cla; 
vx=[0:nx-1]*rep(2)*L(t+1,1)/(nx-1); vy=[0:ny-1]*rep(1)*L(t+1,2)/(ny-1); vz=[0:nz-1]*rep(3)*L(t+1,3)/(nz-1);
[x,y,z]=meshgrid(vx,vy,vz);
rgb=[0 1 0]; rgb2=[.8 1 .8];
meanphi=mean(mean(mean(phi)));
iso1=isosurface(x,y,z,phi,meanphi,phi);
iso2=isocaps(x,y,z,phi,meanphi);
p1 = patch(iso1,'FaceColor',rgb);
p2 = patch(iso2,'FaceColor',rgb2);
n1 = isonormals(x,y,z,phi,p1);
box on; axis equal tight; xlim([0 rep(2)*L(t+1,1)]); ylim([0 rep(1)*L(t+1,2)]); zlim([0 rep(3)*L(t+1,3)]);
view(3);lighting phong; material([.3 .7 .4]); camproj perspective; camlight(-40 ,40);
set(gca,'Color','white'); set(gcf,'Color','white'); xlabel('x'); ylabel('y'); zlabel('z');

%compute area of surface
v1=iso1.vertices;
f1=iso1.faces;
side1=v1( f1(:,1) , :)-v1(f1(:,2) , :);
side2=v1( f1(:,1) , :)-v1(f1(:,3) , :);
normals=cross(side1,side2);
area= 0.5*sum(sum(normals.*normals))
