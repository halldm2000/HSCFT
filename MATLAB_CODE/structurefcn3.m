function structurefcn(fieldname,timeStep)

temp = mgetfieldmpi3d(fieldname,timeStep);
f(:,:)= temp(:,:,1);
[nx,ny]=size(f);
load SCALARS/L.txt;
lx=L(1,1); ly=L(1,2);

%plot morphology
figure(1); clf; surf(f); pbaspect([nx,ny,.25*nx]); 
axis tight; camlight left; lighting phong; shading interp;
title('real space morphology');

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
figure(2); clf; surf(kx,ky,s); pbaspect([nx,ny,.5*nx]); 
axis tight; colormap white; camlight left; lighting phong;
title('scaled scattering function');
xlabel('wavenumber x'); ylabel('wavenumber y'); zlabel('scaled scattering intensity');

figure(3); clf; pcolor(kx,ky,s); pbaspect([nx,ny,.5*nx]); shading flat; colormap(brighten(gray,.5));
xlabel('kx'); ylabel('ky'); title('scattering intensity');

%convert to radial coords
[theta,k,spolar]=cart2pol(kx,ky,s);
figure(4); plot(k,s); axis tight; title('radial slices through scattering data')
kmax = max(max(k));

%bin the data
nbins=nx;
for (i=1:nbins)
    klow = kmax*(i-1)/nbins;
    khigh= kmax*i/nbins;
    ind = find( klow<=k & k<khigh);
    ssum(i) = sum(s(ind));
    area(i) = length(ind);
end
deltaK = kmax/nbins;
kr = ((1:nbins) - 0.5)*deltaK;
sr = ssum./area;

srScaled = sr./( sum(sr)*deltaK );

figure(5); plot(kr,srScaled,'k.-'); title('Scaled Radial Scattering Function vs wavenumber')
figure(6); loglog(kr,srScaled,'k.-'); title('log Scaled Radial Scattering Function vs log wavenumber')

axis tight;


