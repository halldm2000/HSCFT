function m = movie_multisurf8(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    multisurf8(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end