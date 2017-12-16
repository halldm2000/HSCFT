function analyzecontour
[x,y,z]=peaks(200);
figure(1)
[c,h]=contour(x,y,z,[2 2],'k-.')

level = c(1,1)
npoints = c(2,1)
n=2:npoints;
x1=c(1,n+1)
y1=c(2,n+1)

x1=[x1 x1 x1]; y1=[y1 y1 y1];
dx1= x1-circshift(x1,[0 1]); dy1=y1-circshift(y1,[0 1]); ds1= sqrt(dx1.^2 + dy1.^2);
s1=cumsum(ds1); m=max(s1);
step=m*1e-4;
s=0:step:m;
x=interp1(s1,x1,s,'spline');
y=interp1(s1,y1,s,'spline');

dx= x-circshift(x,[0 1]); dy=y-circshift(y,[0 1]); ds= sqrt(dx.^2 + dy.^2); 
dx_c=x-circshift(x,[0 2]); dy_c=y-circshift(y,[0 2]); ds_c=sqrt(dx_c.^2 + dy_c.^2);
dxds = dx./ds; dyds=dy./ds;
d2x = dxds - circshift(dxds,[0 1]); d2y = dyds - circshift(dyds,[0 1]);
d2xds = d2x./ds; d2yds = d2y./ds;

k = (dxds.*d2yds - dyds.*d2xds);
nk = length(k);
k=k(nk/3:2*nk/3)
l=sum(ds)
figure(2); plot(x,y,'r-'); hold on; plot(x1,y1,'.-'); axis equal; hold off;%
figure(3); plot(s1,x1,'b-',s1,y1,'g-'); hold on; plot(s,x,'r-',s,y,'r-'); hold off;
figure(4); plot(s,dxds,s,dyds);
figure(5); plot(s,d2xds,s,d2yds); xlim([m/3,2*m/3]);
figure(6); plot(k)
figure(7); [hy,hx]=hist(k,30); semilogy(hx,hy);

