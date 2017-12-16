function plotscattering(fieldname,timeStep)

temp = mgetfieldmpi3d(fieldname,timeStep);
f(:,:)= temp(:,:,1);
[nx,ny]=size(f);
load SCALARS/L.txt;
lx=L(1,1); ly=L(1,2);

%Compute scaled scattering function
F = fft2(f)/(nx*ny);
S = abs(F).^2;
S(1,1)=0;
S=fftshift(S);
fbar2 = mean(mean(f)).^2
f2bar= mean(mean(f.^2))
s=S/(f2bar - fbar2);
rep=1; s=interp2(s); 
%plot scaled scattering function
[nx,ny]=size(s);
[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(2^rep*lx); ky=fy/(2^rep*ly);

surf(kx,ky,s);  shading flat; colormap(brighten(hot,.3));
view(2); camproj orthographic;
axis equal tight;
xlim([-.6,.6]); ylim([-.6 .6]);