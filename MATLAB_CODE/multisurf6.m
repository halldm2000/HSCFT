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
    if (i==1) hue=zeros(m,n);
    end;
    % choose color
    cfrac = (i-1)/(numLiquids-1);
    ind = find(surfdata>amplitude);
    hue(ind)=2/3*cfrac;

    if (i<=numLiquids) rgb = hsv2rgb([(2/3)*cfrac, 1.0, 1.0]);
    else rgb = hsv2rgb([0 0 0.2 ]);
    end;
    
    % set color to max component
    
    % set amplitude to max concentration
    amplitude=max(amplitude,surfdata);
end

    
[nx,ny]=size(surfdata);
vx=[0:nx-1]*(L(t+1,2)/(nx-1))*rep(1);
vy=[0:ny-1]*(L(t+1,1)/(ny-1))*rep(2);
[x,y]=meshgrid(vx,vy);
    
hsv=zeros(m,n,3);
rgb=zeros(m,n,3);
hsv(:,:,1)=hue;
hsv(:,:,2)=1.0;
hsv(:,:,3)=amplitude;
rgb = hsv2rgb(hsv);
surf(y,x,amplitude,rgb);
shading interp;

%lighting phong;
%material([1 .3 0.05]);
view(2);
%camlight right; camlight left;
%camlight;
pbaspect([rep(2) rep(1) 0.15*rep(1)] )
xlim([0 max(L(:,1))*rep(2)]);
ylim([0 max(L(:,2))*rep(1)]);
zlim([0 1]);
camproj ortho