function m = movie_msurfmpi4(fieldName,startFrame,step,endFrame,rep)
frame=1;
[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});

for i=startFrame:step:endFrame
    figure(1);
    msurfmpi4(fieldName,i,rep);
    axis off;
    time = i*writeInterval;
    zlim([0 4]);
    weight='bold';
        fsize=12;
        fg=[.9 1 .9];
        bg='black'
        text(0,0,4.1,sprintf('t=%.1f',round(time*10)/10),'Color',fg,'FontSize',fsize,'FontWeight',weight,'BackgroundColor',bg);
        drawnow
    m(frame)=getframe;
    frame=frame+1; 
end