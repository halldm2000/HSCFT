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

clear surfdata;
for i=1:numBlocks
    % read the data
    data = mgetfieldmpi3d( sprintf('phi%d',i-1) ,t);
    data2d(:,:)=data(:,:,1);
    phi=repmat(data2d,rep);
    
    %initialize image to white
    if (i==1)
        [m,n]=size(phi);
        red     =ones(m,n);
        green   =ones(m,n);
        blue    =ones(m,n);
    end
    
    % choose color
    cfrac = (i-1)/(numLiquids-1);
    if (i<=numLiquids) rgb = hsv2rgb([(2/3)*cfrac, 1.0, 1.0]);
    else rgb = hsv2rgb([0 0 0.0 ]);
    end;
    
    % override colors
    switch i
        case {1}; rgb=[1 0 0];
        case {2}; rgb=[1 1 0];
        case {3}; rgb=[0 0 0];
    end;
    % apply subtractive color scheme
    invrgb=1-rgb;
    red     = red-phi*invrgb(1);
    green   = green-phi*invrgb(2);
    blue    = blue-phi*invrgb(3);

    % set surface amplitude to max concentration
    amplitude=max(amplitude,phi);
end

[nx,ny]=size(phi);
vx=[0:nx-1]*(L(t+1,2)/(nx-1))*rep(1);
vy=[0:ny-1]*(L(t+1,1)/(ny-1))*rep(2);
[x,y]=meshgrid(vy,vx);
    
clear color;
color(:,:,1)=max(red,0);
color(:,:,2)=max(green,0);
color(:,:,3)=max(blue,0);

% make the surface and set its visual properties
surf(x,y,amplitude-1,color);
lighting phong;
shading interp;
view(2);

pbaspect([rep(2)*ny rep(1)*nx 0.02*rep(1)*nx] )
material dull;
xlim([0 max(L(:,1))*rep(2)]);
ylim([0 max(L(:,2))*rep(1)]);
zlim([0 1]);
camproj ortho
%camlight left; %camlight right;
%camlight left;
camlight headlight;
