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
contourf(vorticity,20)
shading flat;
colormap jet;
axis equal tight;
caxis([-0.01 0.01]);
hold on; contour(solid,1,'k-'); hold off;