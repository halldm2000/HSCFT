function m = movie_msurfmpi(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    mcontmpi(fieldName,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end