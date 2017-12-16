function multisurf(timeStep,rep)

L=load('SCALARS/L.txt');

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
    if (i<=numLiquids) rgb = hsv2rgb([0, 0.0, 0.99-0.80*(i-1)/(numLiquids-1) ]);
    else rgb = hsv2rgb([0 0 0 ]);
    end;
    %surf(surfdata,'FaceColor',rgb,'EdgeColor','black','EdgeAlpha',0.2);
    surf(surfdata,'FaceColor',rgb,'EdgeColor','none');

    hold on;
end
hold off

lighting phong;
material([0.9 0.1 0]);
axis square tight;
view(2);
camlight headlight;
camlight headlight;
pbaspect([rep(2)*L(1,1) rep(1)*L(1,2) 1*L(1,3)]);
camproj ortho;
