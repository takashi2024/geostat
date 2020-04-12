% variog returns ominidirectional residual variogram.
%
% model derives from likfit function.
%
% The semi-variance will not be evaluated (NaN will be returned) if the
% number of lag pairs is less than the value of minPairs.

function [dist_type, variance] = variog(model,dist,X,Y,minPairs,precision)
    dist = round(dist,precision);
    beta = table2array(model.Coefficients(:,1));
    residual = Y - X*beta;
    dist_type = unique(dist);
    Y_dist = squareform(pdist(residual)); 
    N = length(Y);

    variog = zeros(length(dist_type),1);
    lagpairs = zeros(length(dist_type),1);
    for i = 1:length(dist_type)
        target = repmat(dist_type(i),N);
        target = double(dist == target); % the object logical values
        num = sum(sum(target)); % the number of samples lag ij.
        variog(i) = (0.5/num)*sum(sum(((Y_dist).*target).^2)); 
        lagpairs(i) = num;
    end

    variog(variog==0) = NaN;
    variog(lagpairs<minPairs) = NaN;

    variance = variog(:,1);
end