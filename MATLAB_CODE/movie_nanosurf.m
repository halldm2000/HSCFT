function m = movie_nanosurf(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    nanosurf(i,rep); %axis off;
    m(frame)=getframe;
    frame=frame+1; 
end