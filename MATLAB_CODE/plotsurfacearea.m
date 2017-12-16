function plotsurfacearea(field,times,rep)

for i=1:length(times)
    a(i)=surfacearea(field,times(i),rep);
    plot(times(1:i),a(1:i))
end
