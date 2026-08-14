function [xBest,resBest] = DTP2D_multistart_lsqcurvefit(modelFun,x0,xdata,ydata,lb,ub,options,nStarts,use_parallel)

    if nargin < 9
        use_parallel = false;
    end

    rng(1);
    starts = zeros(numel(x0), nStarts+1);
    starts(:,1) = x0(:);
    for k = 1:nStarts
        starts(:,k+1) = lb(:) + rand(numel(x0),1).*(ub(:)-lb(:));
    end

    nTotal = size(starts,2);
    xAll = cell(nTotal,1);
    rAll = inf(nTotal,1);

    if use_parallel
        try
            pool = gcp('nocreate');
            if isempty(pool)
                parpool;
            end

            parfor k = 1:nTotal
                xStart = reshape(starts(:,k),size(x0));
                try
                    [xNow,resNow] = lsqcurvefit(modelFun,xStart,xdata,ydata,lb,ub,options);
                    xAll{k} = xNow;
                    rAll(k) = resNow;
                catch
                    xAll{k} = [];
                    rAll(k) = inf;
                end
            end
        catch ME
            warning('Parallel multi-start failed. Falling back to serial mode. Message: %s', ME.message);
            use_parallel = false;
        end
    end

    if ~use_parallel
        for k = 1:nTotal
            xStart = reshape(starts(:,k),size(x0));
            try
                [xNow,resNow] = lsqcurvefit(modelFun,xStart,xdata,ydata,lb,ub,options);
                xAll{k} = xNow;
                rAll(k) = resNow;
            catch ME
                warning('Start point %d failed: %s', k, ME.message);
                xAll{k} = [];
                rAll(k) = inf;
            end
        end
    end

    [resBest,idxBest] = min(rAll);
    if isinf(resBest) || isempty(xAll{idxBest})
        error('All multi-start fitting attempts failed.');
    end
    xBest = xAll{idxBest};
end
