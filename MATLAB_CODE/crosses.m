function t=crosses(times,data,threshold)
%Find the times at which data crosses the threshold

above = data>threshold;
aboveshift = circshift(above,[0 1]);
i=find(above~=aboveshift);
deltaT = (times(i)-times(i-1))./(data(i)-data(i-1)).*(threshold-data(i-1))
t=times(i)-deltaT