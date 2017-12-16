function m = movie_multisurf(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    frontradius(i,1); %axis off;
    m(frame)=getframe;
    frame=frame+1; 
end