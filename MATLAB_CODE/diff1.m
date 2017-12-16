function diff1(name1,name2,timeStep,rep)

load SCALARS/maxPhi0.txt;
load SCALARS/maxPhi1.txt;

temp = mgetfieldmpi3d(name1,timeStep);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);

temp = mgetfieldmpi3d(name2,timeStep);
temp2d(:,:)= temp(:,:,1);
field2 = repmat(temp2d,rep);
diff=field-field2;
[nx,ny]=size(diff);
if (nx<64) diff=interp2(diff,1); end;

surf(diff);
colormap gray;
shading interp;
axis equal tight;
view(2);
camproj orthographic;
lighting phong;
material([1 0 0]);

max0=max(maxPhi0)
max1=max(maxPhi1)
maxDiff=5*abs(max0-max1)
%caxis([-maxDiff maxDiff])