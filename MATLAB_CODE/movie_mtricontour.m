function m = movie_mtricontour(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    mtricontour(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end