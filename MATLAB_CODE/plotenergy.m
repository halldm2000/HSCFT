function plotenergy(timeStep,rep)

[labels,vals]=textread('STARTUP_FILES/DOMAIN_SIZE_p0.txt','%s\t%s');
cellval  = vals( find(strcmp('numBlocks',labels)) );
numBlocks=str2num(cellval{1})

runTime=load('SCALARS/runTime.txt');
F=load('SCALARS/F.txt');
S=load('SCALARS/S.txt');
Sc0=load('SCALARS/Sc0.txt');
Sc1=load('SCALARS/Sc1.txt');
U=load('SCALARS/U.txt');
U0=load('SCALARS/U0.txt');
U1=load('SCALARS/U1.txt');

figure(3); 
minindex = min( [ length(runTime) length(S) length(F) length(U) ]);

step=1; tsub=runTime(1:step:minindex);
loglog(tsub,  F(1:step:minindex),'-','LineWidth',2,'Color','k'); hold on;
loglog(tsub,  U(1:step:minindex),'-','LineWidth',2,'Color','r'); hold on;
%plot(tsub, -Sc0(1:step:minindex),'--','LineWidth',1,'Color','r'); hold on;
%plot(tsub, -Sc1(1:step:minindex),'-','LineWidth',1,'Color','r'); hold on;
%plot(tsub, -S(1:step:minindex),'-','LineWidth',2,'Color','b'); hold on;

for i=1:numBlocks
    rgb = hsv2rgb([(i-1)/(numBlocks) 1 1]);
    %plot(tsub,-Sc(i,:),'-','LineWidth',1,'Color',rgb);
end

%hold off

xlabel('Time')%,'FontWeight','bold');
ylabel('Energy Components')%,'FontWeight','bold');