function m = movie_msurfmpi5(fieldName,startFrame,step,endFrame,rep)
frame=1;

[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});

for i=startFrame:step:endFrame
    figure(1);
    msurfmpi5(fieldName,i,rep);
    if (i==startFrame) a=getframe(gca);
        d=size(a.cdata)
    end;
    
    zlim([0 4]);
    time = i*writeInterval;
    weight='bold';
    fsize=12;
    fg=[.9 1 .9];
    bg='black'
    text(2,5,3.1,sprintf('t=%.1f',round(time*10)/10),'Color',fg,'FontSize',fsize,'FontWeight',weight,'BackgroundColor',bg);
    
    m(frame)=getframe(gca,[1,1,d(2),d(1)]);
    axis off;
    frame=frame+1; 
end