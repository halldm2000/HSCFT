function o=orientation(fieldname,timeStep)

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
s=interp2(s,0);     
[nx,ny]=size(s);

[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(lx*nx/nx0);        % Compute wavenumbers
ky=fy/(ly*ny/ny0);          

%figure(3); clf; surf(kx,ky,s); view(2); 
% axis square tight; camproj orthographic; shading flat;
% lighting gouraud; camlight; material([1 0 0]); colormap(brighten(gray,.5));
% xlabel('kx'); ylabel('ky'); title('scattering intensity');

[theta,k,spolar]=cart2pol(kx,ky,s); %convert to radial coords
kmax = max(max(k));

%bin the data
deltaTheta = 2*(pi/180);
nbins=(2*pi)/deltaTheta;
thetaLow=-pi:deltaTheta:(pi-deltaTheta);
thetaHigh=(-pi+deltaTheta):deltaTheta:pi;
thetaMid=(thetaLow+thetaHigh)/2;

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

s2theta = s2sum./wedgeArea;
s2thetaMean = mean(s2theta);

a = sqrt(s2thetaMean - sthetaMean.^2);

sx = stheta*(cos(thetaMid).^2)';
sy = stheta*(sin(thetaMid).^2)';
stotal = sum(stheta);
o = (sx-sy)/stotal;

%figure(7); polar(thetaMid,stheta); hold on 
%polar(thetaMid,s2theta,'r'); hold off;
%figure(8); plot(0.5*(thetaHigh+thetaLow),stheta,'k'); pbaspect([2 1 1]);
%xlabel('angle'); ylabel('s(\theta)'); axis tight;

%figure(3); 
%text(0.9*min(min(kx)),0.9*min(min(ky)),max(max(s)),sprintf('a=%0.2e',a),'Color','yellow','FontWeight','bold','FontSize',12);
