function m = movie_msurfmpi5(fieldName1,fieldName2,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    msurfmpi6(fieldName1,fieldName2,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end