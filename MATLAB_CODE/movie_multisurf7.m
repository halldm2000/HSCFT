function m = movie_multisurf7(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    multisurf7(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end