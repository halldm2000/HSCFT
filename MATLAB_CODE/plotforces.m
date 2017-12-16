function multisurf(timeStep,rep)
figure(4);

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1});

%load forces on shearing surfaces
eforce= load('SCALARS/elasticForceMax.txt');
oforce= load('SCALARS/osmoticForceMax.txt');
wforce= load('SCALARS/wallForceMax.txt');
vforce= load('SCALARS/viscousForceMax.txt');

runTime= load('SCALARS/runTime.txt');

minindex = length(runTime);
minindex = min(minindex, length(eforce));
minindex = min(minindex, length(oforce));
minindex = min(minindex, length(wforce));
minindex = min(minindex, length(vforce));

t=runTime(1:minindex);

step=1;
tsub=runTime(1:step:minindex);
eforceSub=eforce(1:step:minindex);
oforceSub=oforce(1:step:minindex);
wforceSub=wforce(1:step:minindex);
vforceSub=vforce(1:step:minindex);

plot(tsub,eforceSub,'r-','LineWidth',2); hold on
plot(tsub,oforceSub,'b-','LineWidth',2); hold on
plot(tsub,wforceSub,'g-','LineWidth',2); hold on
plot(tsub,vforceSub,'k-','LineWidth',2); hold on
hold off;

xlabel('Time')%,'FontWeight','bold');
ylabel('Max Forces')%,'FontWeight','bold');
pbaspect([2 1 1]);
legend('elastic','osmotic','wall','viscous');