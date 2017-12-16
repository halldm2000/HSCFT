function m = movie_multicsurf(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    multicsurf(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end