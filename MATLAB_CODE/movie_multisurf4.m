function m = movie_multisurf4(startFrame,step,endFrame,rep)
frame=1;

[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});

for i=startFrame:step:endFrame
    figure(1);
    multisurf4(i,rep); axis off;
    if (i==startFrame) a=getframe(gca);
        d=size(a.cdata)
    end;
    
    zlim([0 4]);
    time = i*writeInterval;
    weight='bold';
    fsize=12;
    fg=[.9 1 .9];
    bg='black'
    text(1,1,3.1,sprintf('t=%.1f',round(time*10)/10),'Color',fg,'FontSize',fsize,'FontWeight',weight,'BackgroundColor',bg);

    m(frame)=getframe(gca,[0,0,d(2),d(1)]);
    frame=frame+1; 
end