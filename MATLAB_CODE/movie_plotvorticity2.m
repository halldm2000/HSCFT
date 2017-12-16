function m = movie_plotvorticity2(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    plotvorticity2(fieldName,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end