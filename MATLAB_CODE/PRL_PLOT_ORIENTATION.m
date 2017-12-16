cd /Volumes/bifrost.lanl.gov/halldm/DATA/HSCFT_POLYDISPERSE_DATA/Apr212006_132233_TriblockABC111_RectangularPipe_2D64
load orientation; t00=t; o00=o;

cd /Volumes/bifrost.lanl.gov/halldm/DATA/HSCFT_POLYDISPERSE_DATA/Apr212006_163039_TriblockABC111_RectangularPipe_2D64
load orientation; t01=t; o01=o;

cd /Volumes/bifrost.lanl.gov/halldm/DATA/HSCFT_POLYDISPERSE_DATA/Apr242006_144300_TriblockABC111_RectangularPipe_2D64
load orientation; t02=t; o02=o;

cd /Volumes/bifrost.lanl.gov/halldm/DATA/HSCFT_POLYDISPERSE_DATA/Apr242006_145156_TriblockABC111_RectangularPipe_2D64
load orientation; t03=t; o03=o;

cd /Volumes/bifrost.lanl.gov/halldm/DATA/HSCFT_POLYDISPERSE_DATA/Apr242006_144351_TriblockABC111_RectangularPipe_2D64
load orientation; t04=t; o04=o;

cd /Volumes/bifrost.lanl.gov/halldm/DATA/HSCFT_POLYDISPERSE_DATA/Apr242006_144519_TriblockABC111_RectangularPipe_2D64
load orientation; t08=t; o08=o;

figure(4);
semilogx(t00,o00,'s-',t01,o01,'o-',t02,o02,'^-',t03,o03,'>-',t04,o04,'x',t08,o08,'<-');

pbaspect([1.5 1 1]);
xlim([1 90]);
xlabel('Time','FontSize',10); ylabel('Orientation','FontSize',10);
%legend('Wi=0.0048','Wi=0.0125','Wi=0.0285','Wi=0.0416','Wi=0.0532','Wi=0.1074')
legend('Wi=0.005','Wi=0.013','Wi=0.029','Wi=0.042','Wi=0.053','Wi=0.108',0);
legend boxoff
figure(5);
Wi          =[0.0285   ,0.0457  ,0.0532 ,0.0805  ,0.1074 ,0.1314]
timeThresh  =[170.8    ,80.32   ,65.75  ,42.03  ,29.86  ,19]
h=plot(Wi,timeThresh,'ko','MarkerFaceColor','black');
xlabel('Wi','FontSize',18); ylabel('t_{o=40}','FontSize',18);
set(gca,'FontSize',18);
%xlim([0 0.14]);
pbaspect([1.5 1 1])
%box off
%  x0=0.020;
%  x=[0.025+2e-3:1e-3:0.15];
%  alpha=-2.0
%  y=(x-0.005).^(-1.41)
%  plot(x,y,'--');hold off



