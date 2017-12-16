function m = movie_isophi1(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isophi1(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end