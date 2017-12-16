function m = movie_isophis_cutaway(which,startFrame,step,endFrame,rep,cutfracs)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isophis_cutaway(which,i,rep,cutfracs);
    m(frame)=getframe;
    frame=frame+1; 
end