function m = movie_isophis(which,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isophi(which,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end