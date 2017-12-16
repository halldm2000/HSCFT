function msurfmpi7(name1,timeStep1,timeStep0,rep)


temp = mgetfieldmpi3d(name1,timeStep1);
temp2d(:,:)= temp(:,:,1);
field = repmat(temp2d,rep);

temp = mgetfieldmpi3d(name1,timeStep0);
temp2d(:,:)= temp(:,:,1);
field2 = repmat(temp2d,rep);
diff=field-field2;
[nx,ny]=size(diff);

temp = mgetfieldmpi3d('solid',timeStep1);
temp2d(:,:)= temp(:,:,1);
wall = repmat(temp2d,rep);

if (nx<64) diff=interp2(diff,1,'cubic'); end;
if (nx<64) wall=interp2(wall,1,'cubic'); end;
maxdiff=abs(max(max(diff)))
%diff=diff/maxdiff;

surf(diff); hold on;
shading interp;
lighting phong;
material([1 0 0]);
axis tight equal;
view(2);
camproj orthographic;

%caxis([-0.5*maxdiff 0.5*maxdiff])
%caxis([-0.03 0.03])

colormap gray;
%colormap jet;

hold off