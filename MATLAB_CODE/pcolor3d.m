function stereosurf(filename,timeStep,rep)
figure(1); msurfmpi(filename,timeStep,rep); colormap white; material([0.5,0.5,0.5]);

view([-35 30]); set(gcf,'Color','k'); axis off;
f1=getframe(gcf);
figure(2); msurfmpi(filename,timeStep,rep); colormap white; material([0.5,0.5,0.5]);
view([-30 30]); set(gcf,'Color','k'); axis off;
f2=getframe(gcf);
stereo=f2;
stereo.cdata(:,:,1)=f1.cdata(:,:,1);
stereo.cdata(:,:,2)=0;
stereo.cdata(:,:,3)=f2.cdata(:,:,3);
figure(3); image(stereo.cdata); set(gcf,'Color','k');
