function multisurf5(t,rep)
figure(1);
L = load('SCALARS/L.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})
numLiquids = numBlocks-numSolids

amplitude=0;
phase=0;


% read the data
for i=1:numBlocks
    data(i,:,:,:) = mgetfieldmpi3d( sprintf('phi%d',i-1) ,t);
    repdata(i,:,:,:)=repmat(data(i,:,:,:),[1 rep(1) rep(2) 2]);
    phi(i,:,:)=repdata(i,:,:,1);
    surfdata(:,:)=phi(i,:,:);   
    surfdata=surfdata';
    [m,n]=size(surfdata);
    if (i==1)
        red     =zeros(m,n);
        green   =zeros(m,n);
        blue    =zeros(m,n);
    end
    % choose color
    cfrac = (i-1)/(numLiquids-1);
    if (i<=numLiquids) rgb = hsv2rgb([(2/3)*cfrac, 1.0, 1.0]);
    else rgb = hsv2rgb([0 0 0.1 ]);
    end;
    
    % set color to max component
    ind = find(surfdata>amplitude);
    red(ind)=rgb(1);
    green(ind)=rgb(2);
    blue(ind)=rgb(3);
    
    % set amplitude to max concentration
    amplitude=max(amplitude,surfdata);
end

    
[nx,ny]=size(surfdata);
vx=[0:nx-1]*(L(t+1,2)/(nx-1))*rep(1);
vy=[0:ny-1]*(L(t+1,1)/(ny-1))*rep(2);
[x,y]=meshgrid(vx,vy);
    
clear color;
color(:,:,1)=red.*amplitude;
color(:,:,2)=green.*amplitude;
color(:,:,3)=blue.*amplitude;

surf(y,x,amplitude,color);
shading interp;

lighting phong;
material([.7 .3 0.05]);
view(2);
camlight right; camlight left;
camlight;
pbaspect([rep(2) rep(1) 0.10*rep(1)] )
xlim([0 max(L(:,1))*rep(2)]);
ylim([0 max(L(:,2))*rep(1)]);
zlim([0 1]);
camproj ortho