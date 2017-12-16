function R=plotlengthscale(fieldname,timeSteps)

[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1})

for i=1:length(timeSteps)
% load data
temp = mgetfieldmpi3d(fieldname,timeSteps(i));
f(:,:)= temp(:,:,1);
[nx,ny]=size(f);
load SCALARS/L.txt;
lx=L(1,1); ly=L(1,2);

%Compute scaled scattering function
F = fft2(f)/(nx*ny);
S = abs(F).^2;
S(1,1)=0;
S=fftshift(S);
fbar2 = mean(mean(f)).^2;
f2bar= mean(mean(f.^2));
s=S/(f2bar - fbar2);
rep=0;
s=interp2(s,rep);
[nx,ny]=size(s);

[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);
kx = fx/(2^rep*lx); ky=fy/(2^rep*ly);

%convert to radial coords
[theta,k,spolar]=cart2pol(kx,ky,s);
kpeak=sum(sum(k.*s))/sum(sum(s));
R(i)=1./kpeak;
end;
figure(4); plot(timeSteps(1:i)*writeInterval,R.^3,'ks','MarkerSize',5); drawnow;
xlabel('Time','FontSize',12);
ylabel('(1/k_1)^3','FontSize',12);
pbaspect([2 1 1]);
legend off;
