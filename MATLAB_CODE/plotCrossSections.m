function plotCrossSections(t)
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
    phi(i,:,:)=data(i,:,:,1);
    surfdata(:,:)=phi(i,:,:);   
    surfdata=surfdata';
end

[nx,ny]=size(surfdata);
vx=[0:nx-1]*(L(t+1,2)/(nx-1));
vy=[0:ny-1]*(L(t+1,1)/(ny-1));

a=32;
plot(vx,phi(1,:,a),'-k',vx,phi(2,:,a),'--k','LineWidth',1.2)
axis tight; ylim([0 1]);
pbaspect([2 1 1]);
ylabel('\phi');
xlabel('x');
legend('A','B',1); legend boxoff