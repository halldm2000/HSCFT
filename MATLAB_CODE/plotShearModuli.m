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
plot(stress,'k-'); title('stress');
% h = distance between solid walls in the top channel
h = L(1,2)/2;
pbaspect([1 .5 1]);

figure(3); 
% strainRate in the top channel
strain = (rCM0-rCM1)*(1/h+1/h);
plot(strain,'k-'); title('strain');
pbaspect([1 .5 1]);

figure(4);
strainRate = -(vCM0-vCM1)*(1/h+1/h);

% TODO: figure out why I need a factor of 10
viscosity = stress./(strainRate+(strainRate==0))/10;


plot(strain(:,1),stress(:,1)); 

title('Stress vs Strain');
xlabel('strain_{xy}');
ylabel('stress_{xy}');
%axis tight
%ylim([0.3 1]);