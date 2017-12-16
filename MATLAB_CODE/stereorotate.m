function m=stereorotate(filename,timeStep,rep)
frame=1;
for angle=0:5:360
    stereosurf(filename,timeStep,angle,rep);
    m(frame)=getframe;
    frame=frame+1; 
end