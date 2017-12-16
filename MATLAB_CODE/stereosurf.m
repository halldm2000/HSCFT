function stereosurf(filename,timeStep,angle,rep)
figure(1); msurfmpi(filename,timeStep,rep); colormap white; material([0.2, 0.4, 0.2]);

view([-angle-2 20]); set(gcf,'Color','k'); axis off; axis vis3d;
f1=getframe(gcf);
figure(2); msurfmpi(filename,timeStep,rep); colormap white; material([0.2, 0.8 ,0.2]);
view([-angle 20]); set(gcf,'Color','k'); axis off; axis vis3d;
f2=getframe(gcf);
stereo=f2;
stereo.cdata(:,:,1)=f1.cdata(:,:,1);
stereo.cdata(:,:,2)=0;
stereo.cdata(:,:,3)=f2.cdata(:,:,3);
figure(3); image(stereo.cdata); set(gcf,'Color','k');
