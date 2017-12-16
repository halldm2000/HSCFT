function isophis(which,timeStep,rep,cutFractions)
filename='phi0'

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

cellval  = vals( find(strcmp('numSolids',labels)) );
numSolids=str2num(cellval{1});
numLiquids = numBlocks-numSolids;

p=[3 2 1];

%load all the phi scalar fields as a single vector field
for i=1:numLiquids
    temp = mgetfieldmpi3d(sprintf('phi%d',i-1),timeStep);
    temp = interp3(permute(temp,p),rep);
    phi(:,:,:,i) = temp;
    [nz,ny,nx]=size(temp);
    n=size(temp);

end

figure(1); cla;
cutIndex=round(n.*cutFractions);

% plot the isosurfaces for all fields included in the array "which"
for j=1:length(which)

    target = which(j);
    % compute maximum value of all fields not equal to field of interest
    maxOther=zeros(nz,ny,nx);
    for i=1:numLiquids
        if (i~=target) 
            maxOther(:,:,:) = max(phi(:,:,:,i),maxOther(:,:,:));
        end;
    end;

    %clear all data above index cutIndex
    
    diff0(:,:,:) = phi(:,:,:,target)-maxOther(:,:,:);
    diff1(:,:,:) = diff0(1:cutIndex(1),1:cutIndex(2),1:cutIndex(3));
    rgb = hsv2rgb([(2/3)*(target-1)/(numLiquids-1) 1 1 ]);
    p1 = patch(isosurface(diff1,0,diff1),'FaceColor',rgb,'EdgeColor','none');
    p2 = patch(isocaps(diff1,0),'FaceColor',rgb,'EdgeColor','none');
    isonormals(diff1,p1);
end;

%Set shape of the plot
axis tight equal;
%xlim([1 nx]);ylim([1 ny]);zlim([1 nz]);box on;

%Set appearance of the surfaces
material([0.1 1 .4]);
lighting phong;
camproj perspective;
camlight right ;

view([150 30]);
