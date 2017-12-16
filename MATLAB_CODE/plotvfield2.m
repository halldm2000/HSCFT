function plotvfield2(field,timeStep,rep)
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
multisurf4(timeStep,rep); hold on;

step=nx/32;
vx=[0:step:nx-1]*rep(1)*(L(timeStep+1,2)/(nx-1));
vy=[0:step:ny-1]*rep(2)*(L(timeStep+1,1)/(ny-1));
[x,y,z]=meshgrid(vy,vx,3);

% subtract off the mean flow velocity
%fieldx=fieldx-mean(mean(fieldx));
%fieldy=fieldy-mean(mean(fieldy));

% velocity vectors
h=quiver3(x,y,z,fieldx(1:step:nx,1:step:ny),fieldy(1:step:nx,1:step:ny),0*z,1); hold off;
set(h,'Color',[.5 .5 .5],'LineWidth',1.1);
