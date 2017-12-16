load SCALARS/netForce0.txt;
load SCALARS/netForce1.txt;
load SCALARS/vCM0.txt;
load SCALARS/vCM1.txt;
load SCALARS/netTorque1.txt;
load SCALARS/L.txt
load SCALARS/runTime.txt;

stress = (netForce0-netForce1)/L(1,1);
strainRate = -(vCM0-vCM1)/L(1,2);
viscosity = stress./strainRate;
range = 2:length(viscosity);
plot(range/10,viscosity(range,1),'.-')
%semilogy(range/10,viscosity(range,1),'.-')
