function m = movie_msurfmpi3(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    msurfmpi3(fieldName,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end