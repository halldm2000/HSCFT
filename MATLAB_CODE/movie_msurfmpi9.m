function m = movie_msurfmpi(fieldName,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    msurfmpi9(fieldName,i,rep);
    zlim([0 1]);
    m(frame)=getframe;
    frame=frame+1; 
end