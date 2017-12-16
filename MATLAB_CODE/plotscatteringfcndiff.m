function structurefcndiff(fieldname,timeStep)

temp = mgetfieldmpi3d(fieldname,timeStep);
temp2 = mgetfieldmpi3d(fieldname,0);

f(:,:)= temp(:,:,1);%-temp2(:,:,1);
[nx0,ny0]=size(f);
load SCALARS/L.txt;
lx=L(1,1); ly=L(1,2);

%plot morphology
figure(1); clf; surf(f); pbaspect([nx0,ny0,.25*nx0]); 
axis tight; camlight left; lighting phong; shading interp;
title('morphology');

%Compute scaled scattering function
F = fft2(f)/(nx0*ny0);
S = abs(F).^2;
S(1,1)=0;
S=fftshift(S);
fbar2 = mean(mean(f)).^2
f2bar= mean(mean(f.^2))
s=S/(f2bar - fbar2);
s=interp2(s,0); [nx,ny]=size(s);

%plot scaled scattering function
[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(lx*nx/nx0); ky=fy/(ly*ny/ny0);
figure(2); clf; surf(kx,ky,s); pbaspect([nx,ny,.5*nx]); shading flat;
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
    ssum2(i)= sum(s(ind).^2);
    area(i) = length(ind);
end
deltaK = kmax/nbins;
kr = ((1:nbins) - 0.5)*deltaK;
sr = ssum./(area + (area==0));
s2r = ssum2./(area + (area==0));
sr2 = sr.*sr;
sum(sr2)
anisotropy = (s2r - sr2)/(sum(sr2));
srScaled = sr./( sum(sr)*deltaK );
maxindex=find(srScaled==max(srScaled)); krmax=kr(maxindex);

figure(5); plot(kr/krmax,srScaled,'k.-'); title('Scaled Radial Scattering Function vs wavenumber'); axis tight; hold on;
figure(6); loglog(kr/krmax,srScaled,'k.-'); title('log Scaled Radial Scattering Function vs log wavenumber'); axis tight; hold on;
figure(7); semilogx(kr,anisotropy); hold on; title('anisotrpy vs log k');
axis tight;

totalaniostropy = sum(anisotropy*deltaK)
