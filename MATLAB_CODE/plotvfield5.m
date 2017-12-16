function plotvorticity3(field,timeStep,rep)

temp = mgetfieldmpi3d('solid',timeStep);
temp2d(:,:)= temp(:,:,1);
solid = repmat(temp2d,rep);
     
temp = mgetfieldmpi3d(sprintf('%s0',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldx=repmat(temp2d,rep);

temp = mgetfieldmpi3d(sprintf('%s1',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldy=repmat(temp2d,rep);

mag = sqrt(fieldx.^2+fieldy.^2);
vorticity = curl(fieldx,fieldy);
[nx,ny]=size(vorticity);

figure(1);
top=max(max(mag));

surf(mag-top);
view(2); camproj ortho;
shading interp;
lighting phong;
load vortexcmap;
colormap jet;
axis  tight;
pbaspect([size(solid'),1]);
caxis([0 1]-top);

hold on; 
[c,h]=contour(solid,[0.01 0.01],'k-');
[c,h]=contour(solid,[0.5 0.5],'k-');
set(h,'LineWidth',2);
%hold on; h=surf(solid-top-0.01,'FaceColor','black'); %shading interp;

nmin=min(nx,ny)
step=nmin/16;
[x,y,z]=meshgrid(1:step:ny,1:step:nx,0);
h=quiver3(x,y,z,fieldx(1:step:nx,1:step:ny),...
    fieldy(1:step:nx,1:step:ny),0,1); hold off;
set(h,'Color',[0 0 .5]);
hold off;
xlim([1 ny]);
ylim([1 nx]);
zlim([-top 1]);