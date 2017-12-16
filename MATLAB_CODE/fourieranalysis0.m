function fourieranalysis(timeStep)

[kx0,ky0,s0,f0]=fourierSpectra('phi0',timeStep);
[kx1,ky1,s1,f1]=fourierSpectra('phi1',timeStep);
[kx2,ky2,s2,f2]=fourierSpectra('phi2',timeStep);

figure(1); 
contour(f0,1,'r');colormap gray; hold on; 
contour(f1,1,'g'); hold on;
contour(f2,1,'b'); hold off;

figure(3); clf; surf(kx0,ky0,s0); view(2); hold on
axis square tight; camproj orthographic; shading flat;
lighting gouraud; camlight; material([1 0 0]); colormap(brighten(gray,.5));

figure(4);
level1=1e-3;
level2=1e-2;
n=1
[c,h]=contour(kx0,ky0,s0,n,'k'); axis square tight; hold on; set(h,'Color','red');
[c,h]=contour(kx1,ky1,s1,n,'k'); axis square tight; hold on; set(h,'Color','green');
[c,h]=contour(kx2,ky2,s2,n,'k'); axis square tight; hold on; set(h,'Color','blue');

axis square;
set(gca,'Color','black')

xl=xlim; yl=ylim;
figure(5);
set(gcf, 'Renderer','OpenGL');
h1=surf(kx0,ky0,s0); set(h1,'FaceColor','red','EdgeColor','none'); hold on;
h2=surf(kx0,ky0,s1); set(h2,'FaceColor','green','EdgeColor','none'); hold on;
h3=surf(kx0,ky0,s2); set(h3,'FaceColor','blue','EdgeColor','none'); hold on;
base = 1.0e-4*ones(size(s0));
h4=surf(kx0,ky0,base); set(h4,'FaceColor','black','EdgeColor','none'); hold off;
lighting gouraud; axis square;
pbaspect([1 1 1]); 
view(3);
material dull;
camlight; camlight;
axis tight;
xlim([-1 1]); ylim([-1 1]);