function [kx,ky,s,f]=fourierSpectra(fieldname,timeStep);
temp = mgetfieldmpi3d(fieldname,timeStep);
temp2 = mgetfieldmpi3d(fieldname,0);

f(:,:)= temp(:,:,1)-temp2(:,:,1);
[nx0,ny0]=size(f);
load SCALARS/L.txt; lx=L(1,1); ly=L(1,2);

%Compute scaled scattering function
F = fft2(f)/(nx0*ny0);      % Foruier transform
S = abs(F).^2;              % Scattering intensity
S(1,1)=0;                   % Substract off mean scatter
S=fftshift(S);              % Shift 0 to center of plot
fbar2 = mean(mean(f)).^2;   
f2bar= mean(mean(f.^2));
s=S;%/(f2bar - fbar2);        % Compute normalized scattering function
s=interp2(s,2);     
[nx,ny]=size(s);
[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(lx*nx/nx0);        % Compute wavenumbers
ky=fy/(ly*ny/ny0);          