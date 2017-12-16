function plotGinzburgLandau

figure(1);
psi=[-1:0.01:1];
fHot = GL(310,psi);
fCritical= GL(300,psi);
fCold = GL(290,psi);
plot(psi,fHot,'r',psi,fCritical,'g',psi,fCold,'b','LineWidth',1.2)
xlim([-1 1]); ylim([-10 30]);
xlabel('\psi','FontSize',12); ylabel('f(\psi)','FontSize',12);
pbaspect([2 1 1])
text(0.25,20,'T>T_c')
text(0.36,5,'T=T_c')
text(0.52,0,'T<T_c')

figure(2);
q=[0:0.01:2];
a0=0; a2=10; a4=400; Tc=300; 
T=290; Tc=300; psi0=0.5;
f0=a0 + a2*(T-Tc)*psi0^2/2 + a4*psi0^4/4
kappa=1; M=1;
R=-M*f0*(q.^2).*(1 + 2*kappa*(q.^2)/f0);
plot(q,R,'LineWidth',1.2); pbaspect([2 1 1]);
xlabel('wavenumber, q','FontSize',12); ylabel('Amplification factor','FontSize',12); ylim([-5 6]);
i=find(R==max(R));
text(q(i),-4.8,'q_{max}');
line([q(i) q(i)],[-5 6],'LineStyle',':')
%---------------------------
function f=GL(T,psi)
a0=0; a2=10; a4=400; Tc=300; 
f = a0 + a2*(T-Tc)*psi.^2 + a4*psi.^4;