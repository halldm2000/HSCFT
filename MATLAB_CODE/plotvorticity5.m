function plotvorticity5(field,timeStep,rep)

temp = mgetfieldmpi3d('solid',timeStep);
temp2d(:,:)= temp(:,:,1);
solid = repmat(temp2d,rep);
     
temp = mgetfieldmpi3d(sprintf('%s0',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldx=repmat(temp2d,rep);

temp = mgetfieldmpi3d(sprintf('%s1',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldy=repmat(temp2d,rep);

%mag = sqrt(fieldx.^2+fieldy.^2);
vorticity = curl(fieldx,fieldy);
[nx,ny]=size(vorticity);

figure(1);

surf(vorticity);
view(2);
shading interp;
lighting phong;
colormap jet % hot
axis equal tight;
caxis([-0.002 0.002]);
hold on; contour(solid,1,'k-');

step=nx/32;
[x,y]=meshgrid(1:ny,1:nx);
h=streamslice(fieldx,fieldy,2); hold off;
set(h,'Color',[0 0 0]);
% h=quiver3(x,y,z,fieldx(1:step:nx,1:step:ny),...
%     fieldy(1:step:nx,1:step:ny),0*z,1); hold off;
% set(h,'Color','black');
hold off;