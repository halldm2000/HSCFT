function plotshearmoduli

load SCALARS/netForce1.txt;
load SCALARS/netForce2.txt;
load SCALARS/vCM1.txt;
load SCALARS/vCM2.txt;
load SCALARS/rCM1.txt;
load SCALARS/rCM2.txt;

stress = -(netForce1-netForce2);%TODO: divide by "area" of shearing surfaces
strainrate = vCM1-vCM2;
strain = rCM1-rCM2;

% PLOT STRESS VS STRAIN RATE
figure(3); 
plot(strainrate,stress); hold on
xlabel('d\gamma/dt');
ylabel('T_{xy}');
title('stress vs shear strain rate');

maxstress=max(abs(stress(:,1)))
maxstrain=max(strain)

figure(4);
plot(strain/maxstrain,stress(:,1)/maxstress); hold on
title('stress vs shear strain');
xlim([-1.1,1.1]);
ylim([-1.1,1.1]);
axis square;

figure(5); 

plot(strain/maxstrain); hold on
plot(stress(:,1)./maxstress);

% PLOT STRESS AND STRAIN VS TIME, FIND INPHASE AND OUTOFPHASE PARTS
% figure(5); clf;
% maxstress=max(stress(:,1))
% maxstrain=max(strain)
% nstrain = strain/maxstrain;
% min(vCM1(:,1))
% ind = find(min(vCM1(:,1)) == vCM1(:,1))
% subrange=ind(1):ind(2)-1;
% substress = stress(subrange);
% subnstrain = nstrain(subrange);
% plot(substress/maxstress,'k-','LineWidth',2); hold on
% plot(subnstrain,'k--','LineWidth',2);
% 
% product = (substress/maxstress).*subnstrain;
% Gprime = sum(product)/(pi*length(subrange))
% inphase = Gprime*subnstrain;
% outphase = substress/maxstress-inphase;
% plot(inphase,'-r')
% plot(outphase,'-g')
