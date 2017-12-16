function m = movie_isophis0(which,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isophis0(which,i,0,rep);
    m(frame)=getframe;
    frame=frame+1; 
end