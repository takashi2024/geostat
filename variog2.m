% variog2 returns 2-dimensional residual variogram.
%
% model derives from likfit2 function.
%
% The semi-variance will not be evaluated (NaN will be returned) if the
% number of lag pairs is less than the value of minPairs.

function [dist1_type, dist2_type, variance] = variog2(model,dist1,dist2,X,Y,minPairs,precision)
    dist1 = round(dist1,precision);
    dist2 = round(dist2,precision);
    beta = table2array(model.Coefficients(:,1));
    dist1_type = unique(dist1);
    dist2_type = unique(dist2);
    residual = Y - X*beta;
    Y_dist = squareform(pdist(residual)); 
    N = length(Y);

    variance = zeros(length(dist1_type),length(dist2_type));
    for i = 1:length(dist1_type)
        for j = 1:length(dist2_type)
            target1 = repmat(dist1_type(i),N);
            target2 = repmat(dist2_type(j),N);
            target = double(dist1 == target1 & dist2 == target2); % the object logical values
            num = sum(sum(target)); % the number of samples lag ij.
            if num > minPairs
                variance(i,j) = (0.5/num)*sum(sum((Y_dist.*target).^2)); 
            else
                variance(i,j) = NaN; % lag pairs less than minPairs will not be evaluated.
            end

        end
    end
    variance(variance==0) = NaN;
end
