function multisurf(t,rep)
figure(1);
L = load('SCALARS/L.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})

numLiquids = numBlocks-numSolids

for i=1:numBlocks
    data(i,:,:,:) = mgetfieldmpi3d( sprintf('phi%d',i-1) ,t);
    repdata(i,:,:,:)=repmat(data(i,:,:,:),[1 rep(1) rep(2) 2]);
    phi(i,:,:)=repdata(i,:,:,1);
    surfdata(:,:)=phi(i,:,:);   
    
    cfrac = (i-1)/(numLiquids-1);
    if (i<=numLiquids) rgb = hsv2rgb([(2/3)*cfrac, 1.0, 1.0]);
    else rgb = hsv2rgb([0 0 0.2 ]);
    end;
    
    [nx,ny]=size(surfdata);
    vx=[0:nx-1]*(L(t+1,1)/(nx-1));
    vy=[0:ny-1]*(L(t+1,2)/(ny-1));
    [x,y]=meshgrid(vx,vy);
    if (i==1) pcolor(x,y,surfdata); shading interp; colormap gray; end;
    if (i>numLiquids && i<=numBlocks) [c,h]=contour(x,y,surfdata,[0.5 0.5],'-'); set(h,'Color','g','LineWidth',1.0); end;
    hold on;
    caxis([0.2 0.8]);
end
hold off

lighting phong;
material([.5 .5 0.1]);
view(2);
camlight right; camlight headlight;
axis square tight;
xlim([0 max(L(:,1))]);
ylim([0 max(L(:,2))]);
%pbaspect([rep(2)*L(t+1,1) rep(1)*L(t+1,2) .2*L(t+1,1)]);
pbaspect([max(L(:,1)) max(L(:,2)) 2]);
camproj ortho