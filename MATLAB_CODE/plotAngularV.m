function plotangularv(field,timeStep,rep)

temp = mgetfieldmpi3d('solid',timeStep);
temp2d(:,:)= temp(:,:,1);
solid=repmat(temp2d,rep);

temp = mgetfieldmpi3d('vMag',timeStep);
temp2d(:,:)= temp(:,:,1);
vmag=repmat(temp2d,rep,rep);

[nx,ny]=size(vmag);

temp = mgetfieldmpi3d(sprintf('%s0',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldx=repmat(temp2d,rep,rep);

temp = mgetfieldmpi3d(sprintf('%s1',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldy=repmat(temp2d,rep,rep);

mag = sqrt(fieldx.^2+fieldy.^2);
cav=curl(fieldx,fieldy);
maxcav= max(max(cav))
figure(1); clf;

contourf(cav,20); shading flat; hold on;
colormap jet;
%caxis([-1 1]);
contour(solid,3,'b'); hold on;

[X,Y]=meshgrid(1:nx,1:ny);
scale =0.5;
quiver(X,Y,scale*fieldx,scale*fieldy,2,'-k'); 
axis equal tight;
xlim([1 nx]);
ylim([1 ny]);
hold off;