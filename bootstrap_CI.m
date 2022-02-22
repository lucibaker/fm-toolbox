function [CI, avg, std] = bootstrap_CI(data,fcn,varargin)
%computes 95 percent confidence intervals and mean value of fcn(data) using
%     N bootstrapped samples.
%     output is [CI,avg]


noboot=0;
if length(varargin)>0
    if varargin{1}==0
        %don't bootstrap
        noboot=1;
    end
end

N=1000; % number of samples


if noboot %skip bootstrapping
    avg=feval(fcn,data);
    CI(1)=avg;
    CI(2)=avg;
else
    
    if length(data)>2
        xs=bootstrp(N,fcn,data);
        
        xs=sort(xs);
        
        CI(1)=xs(round(0.05*N));
        CI(2)=xs(round(0.95*N));
        avg=mean(xs);
    else
        disp('not enough data for bootstrap')
        avg=[nan];
        CI=[nan nan];
        
    end

end

