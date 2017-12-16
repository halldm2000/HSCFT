function m = movie_isophi0(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isophis0(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end