function m = tsequence_msurfmpi4(field,targetTimes,rep)
%figure(1);
[m,n] = size(targetTimes);
dx1 = 1/(n);
dx2 = 1/(n*1.01);

runTime=load('SCALARS/runTime.txt');
[rlabels,rvals]=textread('STARTUP_FILES/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});

dy = 1/m;
  for r=1:m
    for c=1:n
        index = (r-1)*n +c;
        subplot('Position',[(c-1)*dx1,0.0,dx2,dy]);
        
        step =0;
        time =0;
        
        % find save near deired target time
        step = round(targetTimes(r,c)/writeInterval)
        time = step*writeInterval

        mpcolor(field,step,rep);
        axis off;
        text(0,-2,sprintf('t=%.1f',round(time*10)/10),'Color','black','FontSize',12,'FontWeight','bold');
        drawnow
    end
  end