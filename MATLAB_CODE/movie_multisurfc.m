function m = movie_multisurfg(startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    multisurfc(i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end