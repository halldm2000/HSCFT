function m = movie_structurefcn2(field,startFrame,step,endFrame)
frame=1;
for i=startFrame:step:endFrame
    %figure(1);
    structurefcn2(field,i);
    m(frame)=getframe;
    frame=frame+1; 
end