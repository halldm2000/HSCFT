function a=anisotropy2(fieldname,timeStep)

temp = mgetfieldmpi3d(fieldname,timeStep);
temp2 = mgetfieldmpi3d(fieldname,0);

f(:,:)= temp(:,:,1);
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
s=interp2(s,0);     
[nx,ny]=size(s);

[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(lx*nx/nx0);        % Compute wavenumbers
ky=fy/(ly*ny/ny0);          

figure(3); clf; pcolor(kx,ky,s); pbaspect([nx,ny,.5*nx]); shading flat; colormap(brighten(gray,.5));
xlabel('kx'); ylabel('ky'); title('scattering intensity');

[theta,k,spolar]=cart2pol(kx,ky,s); %convert to radial coords
kmax = max(max(k));

%bin the data
deltaTheta = 2*(pi/180);
nbins=(2*pi)/deltaTheta;
thetaLow=-pi:deltaTheta:(pi-deltaTheta);
thetaHigh=(-pi+deltaTheta):deltaTheta:pi;

clear ssum;
clear s2sum;
for (i=1:nbins)
    ind = find( thetaLow(i)<=theta & theta<=thetaHigh(i));
    ssum(i) = sum(s(ind));
    s2sum(i)= sum(s(ind).^2);
    samplesInWedge=length(ind);
    wedgeArea(i) = samplesInWedge + (samplesInWedge==0);
end

stheta = ssum./wedgeArea;
sthetaMean = mean(stheta);

stheta2 = stheta.^2;
sthetaMean=mean(stheta);
stheta2Mean = mean(stheta2);

a = sqrt(stheta2Mean - sthetaMean^2);
figure(7); polar(thetaHigh,stheta);
figure(8); plot(0.5*(thetaHigh+thetaLow),stheta,'k'); pbaspect([2 1 1]);
xlabel('angle'); ylabel('s(\theta)'); 
