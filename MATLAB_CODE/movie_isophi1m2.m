function m = movie_isophis(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isophi1m2(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end