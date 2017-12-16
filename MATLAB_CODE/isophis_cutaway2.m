function isophis(which,timeStep,smoothVal,rep,cutFracs1,cutFracs2)
filename='phi0'

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1});
numLiquids = numBlocks-numSolids;

L=load('SCALARS/L.txt');

p=[3 2 1];

%load all the phi scalar fields as a single vector field
for i=1:numBlocks
    temp = mgetfieldmpi3d(sprintf('phi%d',i-1),timeStep);
    temp = interp3(permute(temp,p),smoothVal,'cubic');
    temp = repmat(temp,rep);

    phi(:,:,:,i) = temp;
    [nr,nc,nd]=size(temp);
    n=size(temp);
    ny=nr;
    nx=nc;
    nz=nd;
end

figure(1); cla;
cutIndex1=round(n.*cutFracs1);
cutIndex2=round(n.*cutFracs2);

% plot the isosurfaces for all fields included in the array "which"
for j=1:length(which)

    target = which(j);
    % compute maximum value of all fields not equal to field of interest
    maxOther=zeros(ny,nx,nz);
    for i=1:numBlocks
        if (i~=target) 
            maxOther(:,:,:) = max(phi(:,:,:,i),maxOther(:,:,:));
        end;
    end;

    %clear all data above index cutIndex
    
    diff0(:,:,:) = phi(:,:,:,target)-maxOther(:,:,:);
    clear diff1;
    if (j==1) diff1(:,:,:) = diff0(1:cutIndex1(1),1:cutIndex1(2),1:cutIndex1(3));
    else diff1(:,:,:) = diff0(1:cutIndex2(1),1:cutIndex2(2),1:cutIndex2(3));
    end;
    
    if (target<=numLiquids) rgb = hsv2rgb([(2/3)*(target-1)/(numLiquids-1) 1 1 ]);
    %if (target<=numLiquids) rgb = hsv2rgb([1 0 (2/3)*(target-1)/(numLiquids-1)]);

    else rgb = [.2 .2 .2];
    end;
    
    p1 = patch(isosurface(diff1,0,diff1),'FaceColor',rgb,'EdgeColor','none');
    p2 = patch(isocaps(diff1,0),'FaceColor',rgb,'EdgeColor','none');
    isonormals(diff1,p1);
end;

%Set shape of the plot
axis tight; 
xlim([1 nx]); ylim([1 ny]); zlim([1 nz]);box on;
pbaspect([rep(2)*L(3,1) rep(1)*L(1,1) rep(3)*L(2,1)]);
set(gcf,'renderer','opengl');

material([0.2 .8 .1]);
lighting phong;
camproj perspective;
view([-80 20]);
%view([-80 -20]);
camlight headlight;
camlight left ;

%view([150 30]);
