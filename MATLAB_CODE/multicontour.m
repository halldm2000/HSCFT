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

for i=1:numBlocks
    data(i,:,:,:) = mgetfieldmpi3d( sprintf('phi%d',i-1) ,timeStep);
    repdata(i,:,:,:)=repmat(data(i,:,:,:),[1 rep(1) rep(2) 2]);
    phi(i,:,:)=repdata(i,:,:,1);

    surfdata(:,:)=phi(i,:,:);   
    cfrac = (i-1)/(numLiquids-1);
    if (i<=numLiquids) rgb = hsv2rgb([(2/3)*cfrac, 1.0, 1.0]);
    else rgb = hsv2rgb([0 0 0.1 ]);
    end;
    if (i==1) [c,h]=contourf(surfdata,1); hold on; 
        set(h,'EdgeColor',rgb,'LineWidth',2);
    end;
    contour(surfdata,1); hold on;
end
hold off
pbaspect([L0(1) L1(1) L2(1)]);
%lighting phong;
%material([.3 .7 0]);
%axis square tight;
