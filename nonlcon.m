function [c, ceq] = nonlcon(x,cov1,cov2) 

    if strcmp(cov1,'ProdSum')&&strcmp(cov2,'matern')
%         nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        k = x(4);
%         rho1 = x(5);
%         rho2 = x(6);
%         nu1 = x(7);
%         nu2 = x(8);
    elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'exp')
%         nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        k = x(4);
%         rho1 = x(5);
%         rho2 = x(6);
    elseif strcmp(cov1,'ProdSum')&&strcmp(cov2,'sph')
%         nugget = x(1);
        sill1 = x(2);
        sill2 = x(3);
        k = x(4);
%         rho1 = x(5);
%         rho2 = x(6);
    else
        disp('Please set cov1 as SumMetric/ProdSum and cov2 as matern/exp/sph.')
    end
    
    c = k - 1/max([sill1 sill2]);
    ceq = [];
end