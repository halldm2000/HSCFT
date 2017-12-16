load SCALARS/rCM1.txt;
load SCALARS/rCM2.txt;
load SCALARS/rCM3.txt;
load SCALARS/rCM4.txt;
range = length(rCM1);
r=2:range;
plot(rCM1(r,1),rCM1(r,2),'.-',rCM2(r,1),rCM2(r,2),'.-',rCM3(r,1),rCM3(r,2),'.-',rCM4(r,1),rCM4(r,2),'.-');

hold off;
% 	range2=2:length(rCM2)
% 	x2 = rCM2(range2,1); y2 = rCM2(range2,2); z2 = rCM2(range2,3);
% 	x3 = rCM3(range2,1); y3 = rCM3(range2,2); z3 = rCM3(range2,3);
% 	
% 	plot(x2-x2(1),y2-y2(1),...
%          x3-x3(1),y3-y3(1),...
%          x2-x3-(x2(1)-x3(1)),y2-y3-(y2(1)-y3(1)));
