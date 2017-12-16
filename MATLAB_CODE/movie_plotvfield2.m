function m = movie_plotvfield2(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    plotvfield2(fieldName,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end