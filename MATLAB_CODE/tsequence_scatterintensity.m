function m = tsequence_msurfmpi2(field,targetTimes)
figure(1);
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
        step = round(targetTimes(r,c)/writeInterval)
        time = step*writeInterval
        
        plotscatterintensity(field,step);
        axis off; 
        drawnow
    end
  end