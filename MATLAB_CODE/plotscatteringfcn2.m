function [kr,sr]=plotscatteringfcn2(fieldname,timeStep)

% load data
temp = mgetfieldmpi3d(fieldname,timeStep);
f(:,:)= temp(:,:,1);
[nx,ny]=size(f);
load SCALARS/L.txt;
lx=L(1,1); ly=L(1,2);

%plot morphology
%figure(1); clf; surf(f); pbaspect([nx,ny,.25*nx]); 
%axis tight; camlight left; lighting phong; shading interp;
%title('morphology');

%Compute scaled scattering function
F = fft2(f)/(nx*ny);
S = abs(F).^2;
S(1,1)=0;
S=fftshift(S);
fbar2 = mean(mean(f)).^2
f2bar= mean(mean(f.^2))
s=S/(f2bar - fbar2);
rep=1;
s=interp2(s,rep);
[nx,ny]=size(s);


%plot scaled scattering function 3D
[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(2^rep*lx); ky=fy/(2^rep*ly);

% plot scattering intensity 2D
%figure(3); clf; pcolor(kx,ky,s); pbaspect([nx,ny,.25*nx]); shading flat; colormap(brighten(gray,.5));
%xlabel('kx'); ylabel('ky'); title('scattering intensity');

%convert to radial coords
[theta,k,spolar]=cart2pol(kx,ky,s);
%figure(4); plot(k,s); axis tight; title('radial slices through scattering data')
kmax = max(max(k));

%bin the data
nbins=nx;%./2^rep;
deltaK = kmax/nbins;
kr = ((1:nbins) - 0.5)*deltaK;
clear sr;
for (i=1:nbins)
    mask = ones(nx,ny).*exp(-(k-kr(i)).^2/(1*deltaK).^2); 
    areamask=sum(sum(mask));
    mask=mask/areamask;
     sAtK= mask.*s;
     %figure(4); surf(kx,ky,mask);  shading interp; zlim([0 0.01]); drawnow; 
     sr(i) = sum(sum(sAtK));
end

krInterp=min(kr):0.01:max(kr);
srSmoothed = interp1(kr,sr,krInterp,'spline');
figure(5);
plot(krInterp,srSmoothed,'k-'); title('Scaled Radial Scattering Function vs wavenumber');
hold on;
plot(kr,sr,'k.');
axis tight; xlim([0 0.6]); ylim([0 0.025]); 
xlabel('k','FontSize',12); ylabel('s(k,t)','FontSize',12);

kpeak=sum(sr.*kr)/sum(sr)
