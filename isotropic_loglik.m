function [LogLik] = isotropic_loglik(x,dist,X,Y,REML,cov_model) %dist,X,Y

    if strcmp(cov_model,'matern')
        nugget = x(1);
        sill = x(2);
        rho = x(3);
        nu = x(4);
        V = sill * 1/((2^(nu-1))*gamma(nu)) * ((2*sqrt(nu)*dist)/rho).^nu .* besselk(nu,(2*sqrt(nu)*dist)/rho);
        V(find(dist==0)) = sill;
    elseif strcmp(cov_model,'exp')
        nugget = x(1);
        sill = x(2);
        rho = x(3);
        V = sill * exp(-dist/rho);
    else
        disp('Please set cov_model as matern or exp.')
    end
    
    V = V + diag(repelem(nugget, length(Y))); % add nugget in the diag 
    %r = Y - X * inv(X'*inv(V)*X) * (X'*inv(V)*Y);
    r = Y - X / (X'*(V\X)) * (X'*(V\Y));
    [~, p] = size(X);
    if (REML == 1)
        %LogLik = -0.5*(length(Y) - p)*log(2*pi) - 0.5*log(det(V)) - 0.5*log(det(X'*inv(V)*X)) - 0.5*r'*inv(V)*r; %REML
        LogLik = -0.5*(length(Y) - p)*log(2*pi) - 0.5*log(det(V)) - 0.5*log(det(X'*(V\X))) - 0.5*r'*(V\r); %REML
    else
        %LogLik = -0.5*length(Y)*log(2*pi) - 0.5*log(det(V)) - 0.5*r'*inv(V)*r; % ML
        LogLik = -0.5*length(Y)*log(2*pi) - 0.5*log(det(V)) - 0.5*r'*(V\r); % ML
    end
    
    if (LogLik == -Inf) || (LogLik == Inf)
        LogLik = -realmax^0.5;
    end
end