function m = tsequence_multisurf5(targetTimes,rep)
figure(1); clf;
[m,n] = size(targetTimes);
dx1 = 1/(n);
dx2 = 1/(n*1.02);

runTime=load('SCALARS/runTime.txt');
[rlabels,rvals]=textread('STARTUP_FILES/RunParameters.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});

dy = 1/m;
  for r=1:m
    for c=1:n
        index = (r-1)*n +c;
        %subplot('Position',[0,(n-c)*dx1,dy,dx2]);
        subplot('Position',[(c-1)*dx1,0.0,dx2,dy]);
        
        % find save near deired target time
        step = floor(targetTimes(r,c)/writeInterval)
        time = step*writeInterval
        
        multisurf3(step,rep);
        axis off;
        text(0.5, 0.5,1,sprintf('t=%.1f',round(time*10)/10),'Color',[.95 1 .95],...
        'BackgroundColor','black','FontSize',12,'FontWeight','bold');     
        zlim([0 1]);
        drawnow
    end
  end