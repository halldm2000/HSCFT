function m = movie_isophis(field1,field2,startFrame,step,endFrame,rep)
frame=1;
for i=startFrame:step:endFrame
    figure(1);
    isodiff(field1,field2,i,rep);
    m(frame)=getframe;
    frame=frame+1; 
end