function mplotvfield3d(field,timeStep,rep)

temp = mgetfieldmpi3d('solid',timeStep);
solid=repmat(temp,rep);

temp = mgetfieldmpi3d('vMag',timeStep);
vmag=repmat(temp,rep);

[ny,nx,nz]=size(vmag);
 
temp = mgetfieldmpi3d(sprintf('%s0',field),timeStep);
vx=repmat(temp,rep);

temp = mgetfieldmpi3d(sprintf('%s1',field),timeStep);
vy=repmat(temp,rep);

temp = mgetfieldmpi3d(sprintf('%s2',field),timeStep);
vz=repmat(temp,rep);
clf;
set(gcf,'renderer','openGL');
h=streamslice(vx,vy,vz,[nx/2,nx],[],[1,nz/2]);
set(h,'Color',[0 0 0]);

v=sqrt(vx.^2 + vy.^2 + vz.^2);
hold on;
slice(v,[nx/2],[],[nz/2]);
colormap jet;
shading interp;

rgb=[.2 .2 .2];
p1 = patch(isosurface(solid,0.50,solid),'FaceColor',rgb,'EdgeColor','none');
p2 = patch(isocaps(solid,0.50),'FaceColor',rgb,'EdgeColor','none');
material([0.1 1.2 .3]);
isonormals(solid,p1);
lighting phong;

camlight left;
axis equal tight;
xlim([1 nx]);
ylim([1 ny]);
zlim([1 nz]);
box on;
camproj perspective;
view(3)
