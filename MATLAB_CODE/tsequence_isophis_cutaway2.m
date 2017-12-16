function m = tsequence_isophis_cutaway2(which,targetTimes,smoothval,rep,cutFracs1,cutFracs2)
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
        
        % find save near deired target time
        step = floor(targetTimes(r,c)/writeInterval)
        time = step*writeInterval
        isophis_cutaway2(which,step,smoothval,rep,cutFracs1,cutFracs2)
        %text(2,3,sprintf('t=%.1f',round(time*10)/10),'Color','green','FontSize',12,'FontWeight','bold');
        drawnow
    end
  end