function msurfmpi3(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep,rep);

figure(1);
h=surf(field);
view(2);
axis tight;
pbaspect([1 1 0.15]);
lighting phong;
%shading interp;
camlight ;
camproj ortho;%perspective;
colormap white;
material([.2 .8 .3]);
titleString=sprintf('%t=%d',timeStep);
title(titleString);
%view([0 0]);
view([90 0]);

zlim([0 1]);