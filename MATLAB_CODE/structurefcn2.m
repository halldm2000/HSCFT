function qpeak = structurefcn(fieldname,timeStep)
temp = mgetfieldmpi3d(fieldname,timeStep);
field(:,:)= temp(:,:,1);
[nx,ny]=size(field);

fieldbar = mean(mean(field));
fieldbar2=fieldbar*fieldbar;
field2bar= mean(mean(field.*field));

fieldt = fft2(field)/(nx*ny);
Sq = fieldt.*conj(fieldt);
Sq(1,1)=Sq(1,1)-fieldbar2;
sq=Sq/(field2bar-fieldbar2);
sqshift=fftshift(sq);
[fx,fy]=meshgrid(-nx/2:nx/2-1,-ny/2:ny/2-1);

figure(2);
surf(fx,fy,sqshift); shading faceted; view(2);
axis square tight;
shading flat;
colormap gray;

k=sqrt(fx.^2 + fy.^2);
qpeak = sum(sum(k.*sqshift))/sum(sum(sqshift))
