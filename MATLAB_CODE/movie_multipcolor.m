function m = movie_multisurf(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    multipcolor(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end