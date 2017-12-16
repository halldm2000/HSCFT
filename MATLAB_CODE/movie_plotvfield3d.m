function m = movie_plotvfield3d(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    plotvfield3d(fieldName,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end