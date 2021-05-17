function [LogLik] = loglik2(x,dist1,dist2,X,Y,REML,cov1,cov2) 
    
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
    elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'matern')
        nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        k = x(4);
        rho1 = x(5);
        rho2 = x(6);
        nu1 = x(7);
        nu2 = x(8);
        V = sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1)/rho1)...
            + sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2)/rho2)...
            + k * (sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1)/rho1)).*(sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2)/rho2));
        V(dist1==0 & dist2~=0) = sill1...
            + sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2)...
            + k * sill1 .* (sill2 * 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*dist2(dist1==0 & dist2~=0))/rho2));
        V(dist1~=0 & dist2==0) = sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1)...
            + sill2...
            + k * (sill1 * 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*dist1(dist1~=0 & dist2==0))/rho1)).*sill2;
        V(dist1==0 & dist2==0) = sill1 + sill2 + k*sill1*sill2;
    elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'exp')
        nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        k = x(4);
        rho1 = x(5);
        rho2 = x(6);
        V = sill1*exp(-dist1/rho1)+sill2*exp(-dist2/rho2)+k*(sill1*exp(-dist1/rho1)).*(sill2*(exp(-dist2/rho2)));
    elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'sph')
        nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        k = x(4);
        rho1 = x(5);
        rho2 = x(6);
        V = sill1 * (1 - 1.5*dist1/rho1 + 0.5*(dist1/rho1).^3)...
            + sill2 * (1 - 1.5*dist2/rho2 + 0.5*(dist2/rho2).^3)...
            + k * (sill1 * (1 - 1.5*dist1/rho1 + 0.5*(dist1/rho1).^3)).*(sill2 * (1 - 1.5*dist2/rho2 + 0.5*(dist2/rho2).^3));
        V(dist1>rho1 & dist2<=rho2) = sill2 * (1 - 1.5*dist2(dist1>rho1 & dist2<=rho2)/rho2 + 0.5*(dist2(dist1>rho1 & dist2<=rho2)/rho2).^3);
        V(dist1>rho1 & dist2>rho2) = 0;
        V(dist1<=rho1 & dist2>rho2) = sill1 * (1 - 1.5*dist1(dist1<=rho1 & dist2>rho2)/rho1 + 0.5*(dist1(dist1<=rho1 & dist2>rho2)/rho1).^3);
    else
        disp('Please set cov1 as SumMetric/ProdSum and cov2 as matern/exp/sph.')
    end
    
    
    V = V + diag(repelem(nugget, length(Y))); % add nugget in the diag 
    % r = Y - X * inv(X'*inv(V)*X) * (X'*inv(V)*Y);
    r = Y - X * ((X'*(V\X)) \ (X'*(V\Y)));
    [~, p] = size(X);
    if (REML==1)
        % LogLik = -0.5*(length(Y) - p)*log(2*pi) - 0.5*log(det(V)) - 0.5*log(det(X'*inv(V)*X)) - 0.5*r'*inv(V)*r; %REML
        LogLik = -0.5*(length(Y) - p)*log(2*pi) - 0.5*log(det(V)) - 0.5*log(det(X'*(V\X))) - 0.5*r'*(V\r); %REML
    else
        % LogLik = -0.5*length(Y)*log(2*pi) - 0.5*log(det(V)) - 0.5*r'*inv(V)*r; %ML
        LogLik = -0.5*length(Y)*log(2*pi) - 0.5*log(det(V)) - 0.5*r'*(V\r); %ML
    end
    
    if (LogLik == -Inf) || (LogLik == Inf)
        LogLik = -realmax^0.5;
    end
end