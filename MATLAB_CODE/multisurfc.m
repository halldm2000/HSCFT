function multisurf(timeStep,rep)

L0 = load('SCALARS/L0.txt');
L1 = load('SCALARS/L1.txt');
L2 = load('SCALARS/L2.txt');

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1})

numLiquids = numBlocks-numSolids;
rgb=[ .99   .99   .99;
      0.1   0.1   0.1;
      0   0  .99;
      0.99  .99 0;
      0.0 0.99  0;]
  
for i=1:numBlocks
    data(i,:,:,:) = mgetfieldmpi3d( sprintf('phi%d',i-1) ,timeStep);
    repdata(i,:,:,:)=repmat(data(i,:,:,:),[1 rep(1) rep(2) 2]);
    phi(i,:,:)=repdata(i,:,:,1);

    surfdata(:,:)=phi(i,:,:);   
    cfrac = (i-1)/(numLiquids-1);
    
    surf(surfdata,'FaceColor',rgb(i,:),'EdgeColor','none');

    hold on;
end
hold off

lighting phong;
material([.9 .1 0]);
axis square tight;
view(2);
camlight headlight;
camlight headlight;
pbaspect([rep(2)*L0(1) rep(1)*L1(1) 1*L0(1)]);
camproj ortho;
