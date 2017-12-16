function plotvorticity(field,timeStep,rep)
L = load('SCALARS/L.txt');

temp = mgetfieldmpi3d('solid',timeStep);
temp2d(:,:)= temp(:,:,1);
solid = repmat(temp2d,rep);
     
temp = mgetfieldmpi3d(sprintf('%s0',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldx=repmat(temp2d,rep);

temp = mgetfieldmpi3d(sprintf('%s1',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldy=repmat(temp2d,rep);

vorticity = curl(fieldx,fieldy);
v = sqrt(fieldx.^2+fieldy.^2);
[nx,ny]=size(v);

figure(1);
multisurf2(timeStep,rep); hold on;
%pbaspect([ny nx 0.1*nx]);
%view(2);

step=nx/32;
%[x,y,z]=meshgrid(1:step:nx, 1:step:ny, 10)
vx=[0:step:nx-1]*rep(1)*(L(timeStep+1,2)/(nx-1));
vy=[0:step:ny-1]*rep(2)*(L(timeStep+1,1)/(ny-1));
[x,y,z]=meshgrid(vy,vx,3);
h=quiver3(x,y,z,fieldx(1:step:nx,1:step:ny),fieldy(1:step:nx,1:step:ny),0*z,1); hold off;
set(h,'Color','black');