% In directory with A+B phase separating blend...
times=[60,120,240,480,960];
symbol = {'ko','ks','kd','k^','k<','k>','kv'};
mcolor= {'k','w','k','w','k','w'};
clear kr; clear sr;
figure(4); clf;

for i=1:length(times)
    [kri,sri(i,:),kr,sr]=plotscatteringfcn('phi0',times(i));
    figure(4); plot(kr,sr,char(symbol(i)),'MarkerFaceColor',char(mcolor(i)) ); hold on;
    xlim([0 0.5]); ylim([0, 0.06]);
end;
legend('t=60','t=120','t=240','t=480','t=960');
plot(kri,sri,'k-'); hold off; 
xlabel('k','FontSize',12);
ylabel('s(k,t)','FontSize',12);
pbaspect([1.5,1,1]);