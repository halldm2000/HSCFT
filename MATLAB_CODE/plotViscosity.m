load SCALARS/netForce0.txt;
load SCALARS/netForce1.txt;
load SCALARS/vCM0.txt;
load SCALARS/vCM1.txt;
load SCALARS/netTorque1.txt;
load SCALARS/L.txt
load SCALARS/runTime.txt;
load SCALARS/rCM0.txt;
load SCALARS/rCM1.txt;
figure(2);

% stress =Force per unit area on bottom and top walls
Area = 4*L(1,1)*L(1,3)
stress = (netForce0-netForce1)/Area; 

% h = distance between solid walls in the top channel
h = L(1,2)/2;

% strainRate in the top channel
strainRate = -(vCM0-vCM1)*(1/h+1/h);

% TODO: figure out why I need a factor of 10
viscosity = stress./(strainRate+(strainRate==0))/10;

range1 = 2:1:length(viscosity);
range2 = 2:5:length(viscosity);
%plot(range1/10,viscosity(range1,1),'k-'); hold on
%plot(range2/10,viscosity(range2,1),'k.'); hold off
semilogy(range2/10,viscosity(range2,1),'k.')

title('Transient Visocsity vs. Time');
xlabel('time');
ylabel('\eta_s^*');
pbaspect([1 .5 1]);
axis tight
%ylim([0.3 1]);