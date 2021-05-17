% 'likfit2' will fit the 2-dimensional anisotropic model. It returns model
% (the best parameters with the lowest likelihood) and solutions (optimized
% parameter for each run).
%
% Currently, sum-metric (cov1:SumMetric) model with matern, exponential, 
% spherical kernels (cov2:matern, exp, sph) are available. 

function [model, solutions] = likfit2(x0,coords,X,Y,REML,cov1,cov2,Nrun,lower,upper)
    dist1 = squareform(pdist(coords(:,1))); 
    dist2 = squareform(pdist(coords(:,2))); 
    
    opts = optimoptions('fmincon','Algorithm','interior-point');
%     opts = optimoptions('fmincon','Algorithm','sqp');
    ms = MultiStart('Display','iter','FunctionTolerance',1e-6,'PlotFcn',[],...
        'UseParallel',true,'XTolerance',1e-6);
    stpoints = RandomStartPointSet('NumStartPoints',Nrun, ...
        'ArtificialBound',1e4);
    f = @(x)-loglik2(x,dist1,dist2,X,Y,REML,cov1,cov2);
    
    if strcmp(cov1,'ProdSum')
        f2 = @(x)nonlcon(x,cov1,cov2);
    problem = createOptimProblem('fmincon','x0',x0,'objective',f,'options',opts,...
        'lb',lower,...
        'ub',upper,...
        'nonlcon',f2);
    elseif strcmp(cov1,'SumMetric')
        problem = createOptimProblem('fmincon','x0',x0,'objective',f,'options',opts,...
        'lb',lower,...
        'ub',upper);
    end

    tic;
    [x,fval,exitflag,output,solutions] = run(ms,problem,stpoints);
    toc
    
    % Calculate coefficients
    if strcmp(cov1,'SumMetric')&&strcmp(cov2,'matern')
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
    elseif strcmp(cov1,'SumMetric')&&strcmp(cov2,'exp')
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
    elseif strcmp(cov1,'SumMetric')&&strcmp(cov2,'sph')
        nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        sill3 = x(4);
        rho1 = x(5);
        rho2 = x(6);
        rho3 = x(7);
        alpha = x(8);
        dist3 = sqrt(dist1.^2 + (alpha*dist2).^2);
        V = sill1 * (1 - 1.5*dist1/rho1 + 0.5*(dist1/rho1).^3)...
            + sill2 * (1 - 1.5*dist2/rho2 + 0.5*(dist2/rho2).^3)...
            + sill3 * (1 - 1.5*dist3/rho3 + 0.5*(dist3/rho3).^3);
        V(dist1>rho1 & dist2<=rho2 & dist3<=rho3) = sill2 * (1 - 1.5*dist2(dist1>rho1 & dist2<=rho2 & dist3<=rho3)/rho2 + 0.5*(dist2(dist1>rho1 & dist2<=rho2 & dist3<=rho3)/rho2).^3)...
            + sill3 * (1 - 1.5*dist3(dist1>rho1 & dist2<=rho2 & dist3<=rho3)/rho3 + 0.5*(dist3(dist1>rho1 & dist2<=rho2 & dist3<=rho3)/rho3).^3);
        V(dist1>rho1 & dist2>rho2 & dist3<=rho3) = sill3 * (1 - 1.5*dist3(dist1>rho1 & dist2>rho2 & dist3<=rho3)/rho3 + 0.5*(dist3(dist1>rho1 & dist2>rho2 & dist3<=rho3)/rho3).^3);
        V(dist1>rho1 & dist2>rho2 & dist3>rho3) = 0;
        V(dist1<=rho1 & dist2>rho2 & dist3<=rho3) = sill1 * (1 - 1.5*dist1(dist1<=rho1 & dist2>rho2 & dist3<=rho3)/rho1 + 0.5*(dist1(dist1<=rho1 & dist2>rho2 & dist3<=rho3)/rho1).^3)...
            + sill3 * (1 - 1.5*dist3(dist1<=rho1 & dist2>rho2 & dist3<=rho3)/rho3 + 0.5*(dist3(dist1<=rho1 & dist2>rho2 & dist3<=rho3)/rho3).^3);
        V(dist1<=rho1 & dist2>rho2 & dist3>rho3) = sill1 * (1 - 1.5*dist1(dist1<=rho1 & dist2>rho2 & dist3>rho3)/rho1 + 0.5*(dist1(dist1<=rho1 & dist2>rho2 & dist3>rho3)/rho1).^3);
        V(dist1<=rho1 & dist2<=rho2 & dist3>rho3) = sill1 * (1 - 1.5*dist1(dist1<=rho1 & dist2<=rho2 & dist3>rho3)/rho1 + 0.5*(dist1(dist1<=rho1 & dist2<=rho2 & dist3>rho3)/rho1).^3)...
            + sill2 * (1 - 1.5*dist2(dist1<=rho1 & dist2<=rho2 & dist3>rho3)/rho2 + 0.5*(dist2(dist1<=rho1 & dist2<=rho2 & dist3>rho3)/rho2).^3);
%     elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'matern')
%         nugget = x(1);
%         sill1 = x(2);
%         sill2 = x(3);
%         k = x(4);
%         rho1 = x(5);
%         rho2 = x(6);
%         nu1 = x(7);
%         nu2 = x(8);
%         V = sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1)/rho1)...
%             + sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2)/rho2)...
%             + k * (sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1)/rho1)).*(sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2)/rho2));
%         V(dist1==0 & dist2~=0) = sill1...
%             + sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2)...
%             + k * sill1 .* (sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2));
%         V(dist1~=0 & dist2==0) = sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1)...
%             + sill2...
%             + k * (sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1)).*sill2;
%         V(dist1==0 & dist2==0) = sill1 + sill2 + k*sill1*sill2;
%     elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'exp')
%         nugget = x(1);
%         sill1 = x(2);
%         sill2 = x(3);
%         k = x(4);
%         rho1 = x(5);
%         rho2 = x(6);
%         V = sill1*exp(-dist1/rho1)+sill2*exp(-dist2/rho2)+k*(sill1*exp(-dist1/rho1)).*(sill2*(exp(-dist2/rho2)));
%     elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'sph')
%         nugget = x(1);
%         sill1 = x(2);
%         sill2 = x(3);
%         k = x(4);
%         rho1 = x(5);
%         rho2 = x(6);
%         V = sill1 * (1 - 1.5*dist1/rho1 + 0.5*(dist1/rho1).^3)...
%             + sill2 * (1 - 1.5*dist2/rho2 + 0.5*(dist2/rho2).^3)...
%             + k * (sill1 * (1 - 1.5*dist1/rho1 + 0.5*(dist1/rho1).^3)).*(sill2 * (1 - 1.5*dist2/rho2 + 0.5*(dist2/rho2).^3));
%         V(dist1>rho1 & dist2<=rho2) = sill2 * (1 - 1.5*dist2(dist1>rho1 & dist2<=rho2)/rho2 + 0.5*(dist2(dist1>rho1 & dist2<=rho2)/rho2).^3);
%         V(dist1>rho1 & dist2>rho2) = 0;
%         V(dist1<=rho1 & dist2>rho2) = sill1 * (1 - 1.5*dist1(dist1<=rho1 & dist2>rho2)/rho1 + 0.5*(dist1(dist1<=rho1 & dist2>rho2)/rho1).^3);
    else
        disp('Please set cov1 as SumMetric and cov2 as matern/exp/sph. Currently, ProdSum is reserved for product-sum model for cov1')
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
        sigmasq =ssres/(length(Y)-p);
        estimator = 'REML';
    else
        sigmasq = ssres/length(Y);
        estimator = 'ML';
    end
    
    new_sill1 = sill1 * sigmasq; 
    new_sill2 = sill2 * sigmasq;
    
    if strcmp(cov1,'SumMetric')
        new_sill3 = sill3 * sigmasq;
    end
    
    new_nugget = nugget * sigmasq;
    new_C = sigmasq * C;

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
    
    if strcmp(cov1,'SumMetric')&&strcmp(cov2,'matern')
        Result2 = table([new_nugget 0 0]',[new_sill1 new_sill2 new_sill3]',[rho1 rho2 rho3]',[nu1 nu2 nu3]', [NaN NaN alpha]');
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho', 'Nu', 'Alpha'};
        Result2.Properties.RowNames = {'Dim1','Dim2','Joint'};
    elseif strcmp(cov1,'SumMetric')&&strcmp(cov2,'exp')
        Result2 = table([new_nugget 0 0]',[new_sill1 new_sill2 new_sill3]',[rho1 rho2 rho3]', [NaN NaN alpha]');
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho', 'Alpha'};
        Result2.Properties.RowNames = {'Dim1','Dim2','Joint'};
    elseif strcmp(cov1,'SumMetric')&&strcmp(cov2,'sph')
        Result2 = table([new_nugget 0 0]',[new_sill1 new_sill2 new_sill3]',[rho1 rho2 rho3]', [NaN NaN alpha]');
        Result2.Properties.VariableNames = {'Nugget', 'Sill', 'Rho', 'Alpha'};
        Result2.Properties.RowNames = {'Dim1','Dim2','Joint'};
%     elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'matern')
%         Result2 = table([new_nugget 0]',[new_sill1 new_sill2]',[rho1 rho2]', [nu1 nu2]',[k NaN]');
%         Result2.Properties.VariableNames = {'Nugget','Sill','Rho','nu','k'};
%         Result2.Properties.RowNames = {'Dim1','Dim2'};
%     elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'exp')
%         Result2 = table([new_nugget 0]',[new_sill1 new_sill2]',[rho1 rho2]', [k NaN]');
%         Result2.Properties.VariableNames = {'Nugget','Sill','Rho','k'};
%         Result2.Properties.RowNames = {'Dim1','Dim2'};
%     elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'sph')
%         Result2 = table([new_nugget 0]',[new_sill1 new_sill2]',[rho1 rho2]', [k NaN]');
%         Result2.Properties.VariableNames = {'Nugget','Sill','Rho','k'};
%         Result2.Properties.RowNames = {'Dim1','Dim2'};
    end
      
    model = struct;
    model.Description = '2-dimensional spatial LMM';
    model.cov1 = cov1;
    model.cov2 = cov2;
    model.Coefficients = Result1;
    model.GeoVal = Result2;
    model.negLoglik = -fval;
    model.AIC = AIC;
    model.estimator = char(estimator);
    
    if round(sigmasq,1) ~= 1
        disp('Please increase the value of upper bounds for nugget and sill');
    end
end