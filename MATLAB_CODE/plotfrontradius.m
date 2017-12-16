function plotfrontradius(times)
for i=1:length(times)
    r(i)=frontradius(times(i),0);
    plot(times(1:i),r,'k.'); drawnow;
end;