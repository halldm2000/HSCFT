function mplotflowfield2d(field,timeStep,rep)

temp = mgetfieldmpi3d('solid',timeStep);
temp2d(:,:)= temp(:,:,1);
solid=repmat(temp2d,rep);

temp = mgetfieldmpi3d('phi0',timeStep);
temp2d(:,:)= temp(:,:,1);
phi0=repmat(temp2d,rep);

temp = mgetfieldmpi3d('vMag',timeStep);
temp2d(:,:)= temp(:,:,1);
vmag=repmat(temp2d,rep);

[ny,nx]=size(vmag)

temp = mgetfieldmpi3d(sprintf('%s0',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldx=repmat(temp2d,rep);

temp = mgetfieldmpi3d(sprintf('%s1',field),timeStep);
temp2d(:,:)= temp(:,:,1);
fieldy=repmat(temp2d,rep);

mag = sqrt(fieldx.^2+fieldy.^2);

[labels,vals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
cellval  = vals( find(strcmp('L',labels)) );
L=str2num(cellval{1});
[ny,nx]=size(solid)
lengthPerGP = L/(nx-1);
[x,y]=meshgrid(0:nx-1,0:ny-1);

figure(1);

%mcontmpi('solid',timeStep,rep); hold on;
multisurf3(timeStep,rep); hold on;

%fieldx=fieldx-mean(mean(fieldx));
%fieldy=fieldy-mean(mean(fieldy));

h=streamslice(x*lengthPerGP,y*lengthPerGP,fieldx,fieldy,1.5); hold off;
set(h,'Color',[0 0 0]);
axis equal tight;
xlim([0 (nx-1)*lengthPerGP]);
ylim([0 (ny-1)*lengthPerGP]);
hold off;