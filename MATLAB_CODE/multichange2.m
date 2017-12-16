function multichange2(t,rep)
figure(1);
L = load('SCALARS/L.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})

numLiquids = numBlocks-numSolids
clear surfdata
for i=1:numBlocks
    data(i,:,:,:) = mgetfieldmpi3d( sprintf('phi%d',i-1) ,t);
    repdata(i,:,:,:)=repmat(data(i,:,:,:),[1 rep(1) rep(2) 2]);
    phi(i,:,:)=repdata(i,:,:,1);
    surfdata0(:,:)=phi(i,:,:);
    surfdata=interp2(surfdata0,0);
    
    maxsurf = max(max(surfdata));
    minsurf = min(min(surfdata));
    meansurf =mean(mean(surfdata));
    cfrac = (i-1)/(numLiquids-1);
    if (i<=numLiquids) rgb = hsv2rgb([(2/3)*cfrac, 1,maxsurf-minsurf ]);
    else rgb = hsv2rgb([0 0 0.1 ]);
    end;
    
    [nx,ny]=size(surfdata);
    vx=[0:nx-1]*(L(t+1,1)/(nx-1));
    vy=[0:ny-1]*(L(t+1,2)/(ny-1));
    [x,y]=meshgrid(vx,vy);
 	surf(x,y,surfdata-meansurf,'FaceColor',rgb,'EdgeColor','none');
    hold on;
end
hold off
lighting phong;
material([.9 .1 0]);
view(2);
camlight left;% camlight right;
axis square tight;
xlim([0 max(L(:,1))]);
ylim([0 max(L(:,2))]);
%pbaspect([rep(2)*L(t+1,1) rep(1)*L(t+1,2) .2*L(t+1,1)]);
pbaspect([max(L(:,1)) max(L(:,2)) max(L(:,1))]);
camproj ortho