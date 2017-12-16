function m = movie_plotscatter(field,startFrame,step,endFrame)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    plotscatteringfcn2(field,i); hold off;
    m(frame)=getframe;
    frame=frame+1; 
end