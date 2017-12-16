function plotshearmoduli

load SCALARS/netForce0.txt;
load SCALARS/netForce1.txt;
load SCALARS/vCM0.txt;
load SCALARS/vCM1.txt;
load SCALARS/rCM0.txt;
load SCALARS/rCM1.txt;

stress = -(netForce0-netForce1);%TODO: divide by "area" of shearing surfaces
strainrate = vCM0-vCM1;
strain = rCM0-rCM1;

range=1:500;

% PLOT STRESS VS STRAIN RATE
figure(2); clf;
plot(strainrate(range),stress(range));
xlabel('d\gamma/dt');
ylabel('T_{xy}');
title('stress vs shear strain rate');

% PLOT STRESS AND STRAIN VS TIME, FIND INPHASE AND OUTOFPHASE PARTS
figure(3); clf;
maxstress=max(stress(range))
maxstrain=max(strain(range))
nstrain = strain(range)/maxstrain;

ind = find(max(nstrain) == nstrain)
subrange=ind(1):ind(2)-1;
substress = stress(subrange);
subnstrain = nstrain(subrange);
plot(substress/maxstress,'k-','LineWidth',2); hold on
plot(subnstrain,'k--','LineWidth',2);

product = (substress/maxstress).*subnstrain;
Gprime = sum(product)/(pi*length(subrange))
inphase = Gprime*subnstrain;
outphase = substress/maxstress-inphase;
plot(inphase,'-r')
plot(outphase,'-g')
