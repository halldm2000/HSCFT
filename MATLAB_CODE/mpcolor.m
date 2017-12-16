function mpcolor(filename,timeStep,rep)

temp = mgetfieldmpi3d(filename,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep(1),rep(2));
field = interp2(field,2);
runTime=load('SCALARS/runTime.txt');

[ny,nx]=size(field);
h=pcolor(field);
axis equal tight;
colormap gray;
shading flat;
