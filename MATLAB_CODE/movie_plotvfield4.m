function m = movie_plotvfield4(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    plotvfield4(fieldName,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end