% 'likfit' will fit the isotropic (omnidirectional) model. It returns model
% (the best parameters with the lowest likelihood) and solutions (optimized
% parameter for each run).

function [result, solutions] = likfit(x0,coords,X,Y,REML,cov_model,Nrun,lower,upper)
    dist = squareform(pdist(coords));
    
    ms = MultiStart('Display','iter','FunctionTolerance',1e-6,'PlotFcn',[],...
        'UseParallel',true,'XTolerance',1e-6);
    f = @(x)-isotropic_loglik(x,dist,X,Y,REML,cov_model);

    stpoints = RandomStartPointSet('NumStartPoints',Nrun, ...
        'ArtificialBound',1e4);
    problem = createOptimProblem('fmincon','x0',x0,'objective',f,...
        'lb',lower,...
        'ub',upper);

    tic;
    [x,fval,exitflag,output,solutions] = run(ms,problem,stpoints);
    toc
    
    % Calculate coefficients
    if strcmp(cov_model,'matern')
        nugget = x(1);
        sill = x(2);
        rho = x(3);
        nu = x(4);
        V = sill * 1/((2^(nu-1))*gamma(nu)) * ((2*sqrt(nu)*dist)/rho).^nu .* besselk(nu,(2*sqrt(nu)*dist)/rho);
        V(dist==0) = sill;
    elseif strcmp(cov_model,'exp')
        nugget = x(1);
        sill = x(2);
        rho = x(3);
        V = sill * exp(-dist/rho);
    else
        disp('Please set cov_model as matern or exp.')
    end
        V = V + diag(repelem(nugget, length(Y)));
        C = inv(X'*inv(V)*X);
        beta = C*X'*inv(V)*Y;

    % Transformation of Sill and Nugget for variogram from GeoR 
    ivyx = linsolve(V, [Y,X]);
    xivx = ivyx(:,2:end)' * X;
    xivy = ivyx(:,2:end)' * Y;
    yivy = Y' * ivyx(:,1);
    ssres = yivy - 2 * (beta' * xivy) + beta' * xivx * beta; 
    [~,p] = size(X);

    if (REML==1)
        new_sill =ssres/(length(Y)-p);
        estimator = 'REML';
    else
        new_sill = ssres/length(Y);
        estimator = 'ML';
    end

    new_nugget = nugget * new_sill;
    new_C = new_sill * C;

    z_score = beta./sqrt(diag(new_C));
    p_value = 2*normcdf(abs(z_score),'upper'); % two-sided test
    AIC = 2*(p + length(x0) + 1) - 2*-fval; % beta.size + length(x0) + 1

    Result1 = table(beta,sqrt(diag(new_C)),z_score,p_value);
    Result1.Properties.VariableNames = {'Estimate', 'SE', 'Z score','pValue'};
    
    requestIDs = 'X';
    if (p==1)
        Result1.Properties.RowNames = cellstr('(Intercept)'); % Cell Array
    else
        for k = 1 : (p-1)
            requestID{k} = [requestIDs '_' num2str(k,'%d')]; % Cell Array
        end
        Result1.Properties.RowNames = ['(Intercept)', requestID];
    end

    if strcmp(cov_model,'matern')
        Result2 = table(new_nugget,new_sill,rho,nu);
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho','Nu'};
        Result2.Properties.RowNames = {'Omni'};
    elseif strcmp(cov_model,'exp')
        Result2 = table(new_nugget,new_sill,rho);
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho'};
        Result2.Properties.RowNames = {'Omni'};
    end

    result = struct;
    result.Description = 'Iso';
    result.cov_model = cov_model;
    result.Coefficients = Result1;
    result.GeoVal = Result2;
    result.negLoglik = char((strcat({'Negative log-likelihood is '}, num2str(-fval,'%.3f'))));
    result.AIC = char(strcat({'AIC is '}, num2str(AIC,'%.3f')));
    result.estimator = char(strcat({'Estimator is '}, estimator));
end
