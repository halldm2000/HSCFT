function nanosufr2(t,rep)
figure(1);
L = load('SCALARS/L.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})

numLiquids = numBlocks-numSolids

%set(gcf,'Renderer','OpenGl');

data(:,:,:) = mgetfieldmpi3d( 'phi0' ,t);
repdata(:,:,:)=repmat(data(:,:,:),[rep(1) rep(2) 2]);
surfdata(:,:)=repdata(:,:,1);
    
[nx,ny]=size(surfdata);
vx=[0:nx-1]*rep(1)*(L(t+1,2)/(nx-1));
vy=[0:ny-1]*rep(2)*(L(t+1,1)/(ny-1));
[x,y]=meshgrid(vx,vy);

data(:,:,:) = mgetfieldmpi3d('solid' ,t);
repdata(:,:,:)=repmat(data(:,:,:),[rep(1) rep(2) 2]);
solid(:,:)=repdata(:,:,1);

%pcolor(y,x,surfdata');     
h=surf(y,x,surfdata','EdgeColor','none');
shading interp; colormap autumn;
hold on;

h=surf(y,x,4*solid'-2,'EdgeColor','none','FaceColor',[0 0 0]);
caxis([0 1]);
hold off
view(2);
camproj ortho

lighting phong; material([0.7 0.3 0.0]); 
camlight left; camlight headlight;

axis tight equal;

zlim([-1,2]);
