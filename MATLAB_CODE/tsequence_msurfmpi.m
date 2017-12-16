function m = tsequence_msurfmpi5(field,targetTimes,rep)
figure(1);
[m,n] = size(targetTimes);
dx1 = 1/(n);
dx2 = 1/(n*1.01);

runTime=load('SCALARS/runTime.txt');
[rlabels,rvals]=textread('RUN_PARAMETERS/RUN_PARAMETERS.txt','%s\t%s');
writeIntervalCell = rvals(find(strcmp('WRITE_INTERVAL',rlabels)));
appendIntervalCell= rvals(find(strcmp('APPEND_INTERVAL',rlabels)));
writeInterval=str2num(writeIntervalCell{1});
appendInterval=str2num(appendIntervalCell{1});

dy = 1/m;
  for r=1:m
    for c=1:n
        index = (r-1)*n +c;
        subplot('Position',[(c-1)*dx1,0,dx2,dy]);
        
        step =0;
        time =0;
        
        % find save near deired target time
        while(time < targetTimes(r,c))
            step = step + writeInterval;
            time = runTime(step/appendInterval + 1)
        end  
        if (step>0) 
            step=step -writeInterval;
        end
        msurfmpi(field,step,rep);
        axis off;
        drawnow
    end
  end