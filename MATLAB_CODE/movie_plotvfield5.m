function m = movie_plotvfield5(fieldName,startFrame,step,endFrame,rep)
frame=1;
a=getframe;
for i=startFrame:step:endFrame
    figure(1);
    plotvfield5(fieldName,i,rep); 
    axis off;
    if (i==startFrame) a=getframe(gca);
        d=size(a.cdata)
    end;
    m(frame)=getframe(gca,[0,0,d(2),d(1)]);
    frame=frame+1; 
end