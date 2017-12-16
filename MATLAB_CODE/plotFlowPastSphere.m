function plotFlowPastSphere(timeStep)

% generate analytical data
U= 0.01;
a = 0.50;
x0 = 5; y0=10;
[x,y]=meshgrid(0:0.5:10,0:0.5:20);
[theta,r]=cart2pol((x-x0),(y-y0));
magr = 1 - 1.5*(a./r) + 0.5*(a./r).^3;
magt = 1 - 0.75*(a./r) - 0.25*(a./r).^3;
magr=magr.*(r>a);
magt=magt.*(r>a);
ur = -U*sin(theta).*magr;
ut = -U*cos(theta).*magt;
ux = +ur.*cos(theta)- ut.*sin(theta);
uy = +ur.*sin(theta)+ ut.*cos(theta);

% load numerical results
temp = mgetfieldmpi3d('v0',timeStep);
fieldx(:,:)= temp(:,:,1);

temp = mgetfieldmpi3d('v1',timeStep);
fieldy(:,:)= temp(:,:,1);

[labels,vals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
cellval  = vals( find(strcmp('L',labels)) );
L=str2num(cellval{1});
[ny,nx]=size(fieldx)
lengthPerGP = L/(nx-1);
[x2,y2]=meshgrid(0:nx-1,0:ny-1);
x2=x2.*lengthPerGP;
y2=y2.*lengthPerGP;
vx = interp2(x2,y2,fieldx,x,y);
vy = interp2(x2,y2,fieldy,x,y);
vx=vx.*(r>a);
vy=vy.*(r>a);

% plot them for comparison
clf;
%quiver(x,y,ux,uy,'b'); hold on 
%quiver(x,y,vx,vy,'r'); 
h1=streamslice(x,y,ux,uy); hold on
h2=streamslice(x,y,vx,vy); hold off;
set(h1,'Color','black');
set(h2,'Color','red');
axis equal tight;
