function plotvorticity(field,timeStep,rep)

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
lighting phong;
shading interp;
colormap jet;
axis tight;
pbaspect([ny nx 0.1*nx]);
view(2);
camlight;
material([1 0 0])
max(max(vorticity))
%caxis([-0.002 0.002]);
hold on; [c,h]=contour3(solid,[0.5 0.5],'k-'); hold off;
set(h,'LineWidth',2);
zlim([0 1]);