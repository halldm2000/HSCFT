function m = movie_nanosurf3(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    nanosurf3(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end