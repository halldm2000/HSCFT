load SCALARS/aspectRatioXY.txt;
load SCALARS/aspectRatioXZ.txt;
load SCALARS/volume.txt;
load SCALARS/runTime.txt
figure(6);
plot(runTime,aspectRatioXY-aspectRatioXY(1),'k-'); hold on;
plot(runTime,aspectRatioXZ-aspectRatioXZ(1),'k--'); hold on;
%plot(runTime,aspectRatioXY,'k-'); hold on;
%plot(runTime,aspectRatioXZ,'k--'); hold on;
plot(runTime,volume-volume(1),'b-'); hold on;

load SCALARS/dFdAxy.txt;
load SCALARS/dFdAxz.txt;
load SCALARS/dFdV.txt;
plot(runTime,dFdAxy,'r-'); hold on;
plot(runTime,dFdAxz,'r--'); hold on;
plot(runTime,dFdV,'g-'); hold on;
