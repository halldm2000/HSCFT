function m = movie_multisurf5(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    multisurf5(i,rep); 
    axis off;
    material([.4 .6 .2]);
    m(frame)=getframe;
    frame=frame+1; 
end