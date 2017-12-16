function m = movie_isophis_cutaway(which,startFrame,step,endFrame,smoothVal,rep,cutfracs1,cutfracs2)
frame=1;
for i=startFrame:step:endFrame
    i
    figure(1);
    isophis_cutaway2(which,i,smoothVal,rep,cutfracs1,cutfracs2);
    m(frame)=getframe;
    frame=frame+1; 
end