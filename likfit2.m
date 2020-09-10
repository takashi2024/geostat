% 'likfit2' will fit the 2-dimensional anisotropic model. It returns model
% (the best parameters with the lowest likelihood) and solutions (optimized
% parameter for each run).
%
% Currently, sum-metric model with matern or exponential kernel is available. 

function [model] = likfit2(x0,coords,X,Y,REML,cov_model,Nrun,lower,upper)
    dist1 = squareform(pdist(coords(:,1))); 
    dist2 = squareform(pdist(coords(:,2))); 
    
    ms = MultiStart('Display','iter','FunctionTolerance',1e-6,'PlotFcn',[],...
        'UseParallel',true,'XTolerance',1e-6);
    f = @(x)-sum_metric_loglik(x,dist1,dist2,X,Y,REML,cov_model);

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
        sill1 = x(2);
        sill2 = x(3);
        sill3 = x(4);
        rho1 = x(5);
        rho2 = x(6);
        rho3 = x(7);
        nu1 = x(8);
        nu2 = x(9);
        nu3 = x(10);
        alpha = x(11);
        dist3 = sqrt(dist1.^2 + (alpha*dist2).^2);
        V = sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1)/rho1)...
            + sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2)/rho2)...
            + sill3 * 1/((2^(nu3-1))*gamma(nu3)) * ((2*sqrt(nu3)*dist3)/rho3).^nu3 .* besselk(nu3,(2*sqrt(nu3)*dist3)/rho3);
        V(dist1==0 & dist2~=0) = sill1 + sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2)...
            + sill3 * 1/((2^(nu3-1))*gamma(nu3)) * ((2*sqrt(nu3)*dist3(dist1==0 & dist2~=0))/rho3).^nu3 .* besselk(nu3,(2*sqrt(nu3)*dist3(dist1==0 & dist2~=0))/rho3);
        V(dist1~=0 & dist2==0) = sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1)...
            + sill2 + sill3 * 1/((2^(nu3-1))*gamma(nu3)) * ((2*sqrt(nu3)*dist3(dist1~=0 & dist2==0))/rho3).^nu3 .* besselk(nu3,(2*sqrt(nu3)*dist3(dist1~=0 & dist2==0))/rho3);
        V(dist1==0 & dist2==0) = sill1 + sill2 + sill3;
    elseif strcmp(cov_model,'exp')
        nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        sill3 = x(4);
        rho1 = x(5);
        rho2 = x(6);
        rho3 = x(7);
        alpha = x(8);
        dist3 = sqrt(dist1.^2 + (alpha*dist2).^2);
        V = sill1 * exp(-dist1/rho1) + sill2 * exp(-dist2/rho2) + sill3 * exp(-dist3/rho3);
    else
        disp('Please set cov_model as matern or exp.')
    end
    V = V + diag(repelem(nugget, length(Y)));
    % C = inv(X'*inv(V)*X);
    C = inv(X'*(V\X));
    beta = C*X'*(V\Y);

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
    
    new_sill1 = sill1 * new_sill; %/(sill1 + sill2 + sill3);
    new_sill2 = sill2 * new_sill; %/(sill1 + sill2 + sill3);
    new_sill3 = sill3 * new_sill; %/(sill1 + sill2 + sill3);
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
        Result2 = table([new_nugget NaN NaN]',[new_sill1 new_sill2 new_sill3]',[rho1 rho2 rho3]',[nu1 nu2 nu3]', [NaN NaN alpha]');
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho', 'Nu', 'Alpha'};
        Result2.Properties.RowNames = {'Dim1','Dim2','Joint'};
    elseif strcmp(cov_model,'exp')
        Result2 = table([new_nugget NaN NaN]',[new_sill1 new_sill2 new_sill3]',[rho1 rho2 rho3]', [NaN NaN alpha]');
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho', 'Alpha'};
        Result2.Properties.RowNames = {'Dim1','Dim2','Joint'};
    end
      
    model = struct;
    model.Description = 'Aniso';
    model.cov_model = cov_model;
    model.Coefficients = Result1;
    model.GeoVal = Result2;
    model.negLoglik = char((strcat({'Negative log-likelihood is '}, num2str(-fval,'%.3f'))));
    model.AIC = char(strcat({'AIC is '}, num2str(AIC,'%.3f')));
    model.estimator = char(strcat({'Estimator is '}, estimator));
end
