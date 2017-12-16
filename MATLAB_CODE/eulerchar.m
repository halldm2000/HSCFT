function normalalizedArea =surfacearea(field,t,rep)

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

boxarea = 2*( L(t+1,1)*L(t+1,2) + L(t+1,2)*L(t+1,3) + L(t+1,1)*L(t+1,3) )
%compute area of surface (but not the caps, since domain is periodic)
v1=iso1.vertices;
f1=iso1.faces;
e1=v1( f1(:,1) , :)-v1(f1(:,2) , :);
e2=v1( f1(:,1) , :)-v1(f1(:,3) , :);
e3=v1( f1(:,2) , :)-v1(f1(:,3) , :);
e=[e1; e2; e3];
%normals=cross(side1,side2);
%area= 0.5*sum(sum(normals.*normals))
%normalalizedArea = area/boxarea

vsorted=sort(v1')';
vunique=unique(vsorted,'rows');

esorted=sort(e')';
eunique=unique(e,'rows');
funique=unique(f1,'rows');

numv=length(vunique)
numf=length(funique)
nume=length(eunique)
echar = numf+ numv- nume