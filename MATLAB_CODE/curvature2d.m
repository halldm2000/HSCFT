function curvature2d(filename,t);
data1 = mgetfieldmpi3d(filename,t);
data0 = mgetfieldmpi3d(filename,0);
data=data1-data0;
f=mean(mean(data));

L = load('SCALARS/L.txt');
[nx,ny]=size(data);
vx=[0:nx-1]*(L(t+1,1)/(nx-1));
vy=[0:ny-1]*(L(t+1,2)/(ny-1));
[x,y]=meshgrid(vx,vy);
    
figure(1); 
contourf(x,y,data,1); hold on; colormap gray; axis equal tight;
[c,h]=contour(x,y,data,1,'r-.'); hold off;

i=1; 
ktotal=[]; ltotal=0;
while(i<length(c))
    level = c(1,i)
    npoints = c(2,i)
    n=(i+2):(i+npoints);
    x1=c(1,n);
    y1=c(2,n);
    [k,l]=analyzecontour(x1,y1);
    ktotal=[ktotal k]; ltotal=ltotal+l;
    i=i+npoints+1;
    nk = length(ktotal)
    pause
end;
ltotal
kl=ktotal*ltotal;
%figure(9); [hy,hx]=hist(kl,[-30:0.5:30]);plot(hx,hy/nk,'k'); xlim([-30 30]); hold on;
figure(9); semilogy(kl);