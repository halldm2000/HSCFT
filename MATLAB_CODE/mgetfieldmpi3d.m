function datatrans=mgetfieldmpi3d(fieldName,timeStep)

filename = sprintf('%s/%s_p%d_%d.txt',fieldName,fieldName,0,timeStep);
rawData = importdata(filename);

domainFile= sprintf('STARTUP_FILES/DOMAIN_SIZE_p%d.txt',0);
[dlabels,dvals]=textread(domainFile,'%s\t%f');
nx  = dvals( find(strcmp('nxLocal',dlabels)) );
ny  = dvals( find(strcmp('nyLocal',dlabels)) );
nz  = dvals( find(strcmp('nzLocal',dlabels)) );
numProcesses = dvals( find(strcmp('numProcesses',dlabels)) );

%data stored in format where x changes most slowly
data    = reshape(rawData, nz,ny,nx); 

for i=1:numProcesses-1
    filename = sprintf('%s/%s_p%d_%d.txt',fieldName,fieldName,i,timeStep);
    rawData = importdata(filename);

    domainFile= sprintf('STARTUP_FILES/DOMAIN_SIZE_p%d.txt',i);
    [dlabels,dvals]=textread(domainFile,'%s\t%f');
    nxl  = dvals( find(strcmp('nxLocal',dlabels)) );
    nyl  = dvals( find(strcmp('nyLocal',dlabels)) );
    nzl  = dvals( find(strcmp('nzLocal',dlabels)) );
    nx=nx+nxl;
    
    block = reshape(rawData, nzl,nyl,nxl); 
    [dz,dy,dx]=size(data);
    [bz,by,bx]=size(block);
    data(:,:,dx+1:dx+bx)=block;
end

% coords (x,y,z) corresponds to matrix indices ('column','row','depth')
for x=1:nx
    for y=1:ny
        for z=1:nz
            datatrans(y,x,z)=data(z,y,x);
        end
    end
end