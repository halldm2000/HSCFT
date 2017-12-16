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


for i=1:numBlocks
    % read the data
    data(i,:,:,:) = mgetfieldmpi3d( sprintf('phi%d',i-1) ,t);
    repdata(i,:,:,:)=repmat(data(i,:,:,:),[1 rep(1) rep(2) 2]);
    phi(i,:,:)=repdata(i,:,:,1);
    surfdata(:,:)=phi(i,:,:);   
    surfdata=surfdata';
    [m,n]=size(surfdata);
    
    %initialize image to white
    if (i==1)
        red     =ones(m,n);
        green   =ones(m,n);
        blue    =ones(m,n);
    end
    
    % choose color
    cfrac = (i-1)/(numLiquids-1);
    if (i<=numLiquids) rgb = hsv2rgb([1, 0.0, 1-1.2*cfrac]);
    else rgb = hsv2rgb([0 0 0 ]);
    end;

    if (i==1) rgb=hsv2rgb([1 0 1]); end;
    if (i==2) rgb=hsv2rgb([1 0 .8]); end;
    if (i==3) rgb=hsv2rgb([1 0 0]); end;
        
    
    % apply subtractive color scheme
    invrgb=1-rgb;
    red     = red-surfdata*invrgb(1);
    green   = green-surfdata*invrgb(2);
    blue    = blue-surfdata*invrgb(3);

    % set surface amplitude to max concentration
    amplitude=max(amplitude,surfdata);
end

[nx,ny]=size(surfdata);
vx=[0:nx-1]*(L(t+1,2)/(nx-1))*rep(1);
vy=[0:ny-1]*(L(t+1,1)/(ny-1))*rep(2);
[x,y]=meshgrid(vx,vy);
    
clear color;
color(:,:,1)=max(red,0);
color(:,:,2)=max(green,0);
color(:,:,3)=max(blue,0);

% make the surface and set its visual properties
surf(y,x,amplitude,color);
view(2);
lighting phong;
shading interp;

%material([0.4 .6 0.05])
%pbaspect([rep(2) rep(1) 0.05*rep(1)] )
%camlight(1,0,'infinite');
%camlight left;
%camlight right;
%camlight;
%camlight headlight;
pbaspect([rep(2) rep(1) 0.05*rep(1)] )
material([0.8 0.2 0.0]);
%camlight;

xlim([0 max(L(:,1))*rep(2)]);
ylim([0 max(L(:,2))*rep(1)]);
zlim([0 1]);
camproj ortho