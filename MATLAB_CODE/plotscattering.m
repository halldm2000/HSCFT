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

%plot scaled scattering function
[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/lx; ky=fy/ly;

figure(3); clf; pcolor(kx,ky,s); pbaspect([nx,ny,.5*nx]); shading flat; colormap(brighten(gray,.5));
xlabel('kx'); ylabel('ky'); title('scattering intensity');


