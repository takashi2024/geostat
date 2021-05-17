%% Data import 
% Three dataset (Field1–3) are available for this demo. 
 data = readtable('Projects/GitHub/Data_for_2.5m_mesh/Field1.csv');
% data = readtable('Projects/GitHub/Data_for_2.5m_mesh/Field2.csv');
% data = readtable('Projects/GitHub/Data_for_2.5m_mesh/Field3.csv');

%% Data preparation

% Fixed design matrix (Independent variables)
X = horzcat(repelem(1,length(data.Yield))', data.D1); % The 1st column is for intecept. The 2nd column is for D1–4. 

% Response variable
Y = data.Yield; 

% coordinates
coords = [data.x data.y];

%% Yield and treatment map
figure;
subplot(1,2,1);
scatter(data.x,data.y,30,Y,'filled')
title('Yield map')
xlabel('x coords (m)')
ylabel('y coords (m)')
colormap(gca,parula(64));
cb = colorbar('northoutside');
cb.Label.String = 't ha^{-1}';
pbaspect([1 2 1])

subplot(1,2,2);
scatter(data.x,data.y,30,X(:,2),'filled')
title('Treatment map')
xlabel('x coords (m)')
ylabel('y coords (m)')
colormap(gca,parula(2));
cb = colorbar('northoutside');
pbaspect([1 2 1])

%% Hypothetical treatment effect (0.3 t/ha) is added by Gaussian random number generator
% If the objective of this demo is to estimate the simulated Type I rates,
% it is not necessary to run this section.
rng('default') % For reproducibility
effect = normrnd(0.3,0.1,[sum(X(:,2),1),1]);
meanval = mean(effect);
fprintf('The hypothetical treatment effect is %.4f t/ha \n', meanval)
Y(X(:,2)==1) = Y(X(:,2)==1) + effect; 

%% Fit ordinary least squares (OLS) regression model
result_ols = fitlm(X(:,2:end), Y);

%% Fit isotropic model
rng default % For reproducibility
x0 = [0.5 0.5 10]; % Initial values for parameters (nugget, sill, and rho).
Nrun = 10; % Iterations for optimazation
lower = [1e-9 1e-9 1e-9]; % Lower bound list for the parameters
upper = [1 4 50];% Upper bound list for the parameters
[model_1] = likfit(x0,coords,X,Y,1,'exp', Nrun,lower,upper); % Fit model

%% Fit Anisotropic model
rng default % For reproducibility
x0 = [0.5 0.5 0.5 0.5 1 1 1 1]; % Initial values for parameters (nugget, sill1, sill2, sill3, rho1, rho2, rho3, and alpha)
Nrun = 10; % Iterations for optimazation
lower = [1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9 1e-9]; % Lower bound list for the parameters 
upper = [1 4 4 4 50 100 50 1e2];  % Upper bound list for the parameters 
[model_2] = likfit2(x0,coords,X,Y,1,'SumMetric','exp',Nrun,lower,upper); % Fit model

%% Summary: All results (OLS, isotropic and anisotropic models) are shown.
% OLS model
disp('OLS:')
disp(result_ols)
fprintf('pValue as Z score = %.6f \n\n', ...
    2*normcdf(abs(table2array(result_ols.Coefficients(2,3))),'upper'))

% Isotropic model
names = fieldnames(model_1);
for i = 1:length(names)
    iField = names{i};
    disp(model_1.(iField))
end
fprintf("\n")

% Anisotropic model
names = fieldnames(model_2);
for i = 1:length(names)
    iField = names{i};
    disp(model_2.(iField))
end

%% Plot (residual) experimental and fitted variogram (omnidirectional)

[dist_type, variance_1] = variog(model_1,coords,X,Y,30,0);

% Plot experimental variogram
figure('units','inch','position',[0 0 4.5 3.5]);
plot(dist_type,variance_1,'.','MarkerSize',12);
hold on

% Plot fitted variogram
h = 1:0.1:max(max(dist_type));
nugget = table2array(model_1.GeoVal(:,1));
sill = table2array(model_1.GeoVal(:,2));
rho = table2array(model_1.GeoVal(:,3));

if strcmp(model_1.cov_model,'exp')
    r = nugget + sill*(1-exp(-h/rho)); % For Exponential
elseif strcmp(model_1.cov_model,'sph')
    r = nugget + sill*(1.5*h/rho - 0.5*(h/rho).^3); % For Spherical
    r(h>rho) = nugget + sill; % For Spherical
elseif strcmp(model_1.cov_model,'matern')
    nu = table2array(model_1.GeoVal(:,4)); % For Matern
    r = nugget + sill * (1 - 1/((2^(nu-1))*gamma(nu)) * ((2*sqrt(nu)*h)/rho).^nu .* besselk(nu,(2*sqrt(nu)*h)/rho)); % For matern
end
plot(h,r,'Color','R','linewidth',2)
hold off
title('Variogram')
xlabel('Lag (m)')
ylabel('Semi-variance (t ha^{-1})^{2}')
set(gca,'fontsize', 18,  'FontName', 'Times New Roman','linewidth',1.5);

%% Plot (residual) experimental and fitted variogram (2 dimensional)
[dist1_type,dist2_type,variance_2] = variog2(model_2,coords,X,Y,30,1);

% Plot experimental variogram
figure;
mesh(dist1_type,dist2_type,variance_2','FaceLighting','gouraud','LineWidth',0.5)
title('Experimental variogram')
xlabel('{\it x} lag (m)')
ylabel('{\it y} lag (m)')
zlabel('Semi-variance (t ha^{-1})^{2}')

%% Plot fitted variogram (2 dimensional)
dist1_type_list = 0:2.5:round(max(coords(:,1)));
dist2_type_list = 0:2.5:round(max(coords(:,2)));

m = length(dist1_type_list);
n = length(dist2_type_list);

mat_variog = zeros(m*n,3);

mat_variog(:,1) = repmat(dist1_type_list, 1, n);

ini = 1;
for i = 1:n
    mat_variog(ini:ini+m-1,2) = repmat(dist2_type_list(i),1,m);
    ini = ini + m;
end

%mat_variog(mat_variog==0) = NaN;

xlag = mat_variog(:,1);
ylag = mat_variog(:,2);

[Xlag,Ylag] = meshgrid(xlag, ylag);
Xlag(Xlag==max(xlag)) = NaN; % avoid to connect the initial and last points

nugget = table2array(model_2.GeoVal(1,1));
sill1 = table2array(model_2.GeoVal(1,2));
sill2 = table2array(model_2.GeoVal(2,2));
rho1 = table2array(model_2.GeoVal(1,3));
rho2 = table2array(model_2.GeoVal(2,3));

if strcmp(model_2.cov1,'SumMetric')&&strcmp(model_2.cov2,'matern')
    sill3 = table2array(model_2.GeoVal(3,2));
    rho3 = table2array(model_2.GeoVal(3,3));
    nu1 = table2array(model_2.GeoVal(1,4));
    nu2 = table2array(model_2.GeoVal(2,4));
    nu3 = table2array(model_2.GeoVal(3,4));
    alpha = table2array(model_2.GeoVal(3,5));
    Zlag = sqrt(Xlag.^2 + (alpha*Ylag).^2);
    matern1 = 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*Xlag)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*Xlag)/rho1);
    matern2 = 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*Ylag)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*Ylag)/rho2);
    matern3 = 1/((2^(nu3-1))*gamma(nu3)) * ((2*sqrt(nu3)*Zlag)/rho3).^nu3 .* besselk(nu3,(2*sqrt(nu3)*Zlag)/rho3);
    fitVariog = nugget+sill1+sill2+sill3-(sill1*matern1+sill2*matern2+sill3*matern3);
elseif strcmp(model_2.cov1,'SumMetric')&&strcmp(model_2.cov2,'exp')
    sill3 = table2array(model_2.GeoVal(3,2));
    rho3 = table2array(model_2.GeoVal(3,3));
    alpha = table2array(model_2.GeoVal(3,4));
    fitVariog = nugget + sill1 + sill2 + sill3...
        - (sill1 * exp(-Xlag/rho1)...
        + sill2 * exp(-Ylag/rho2)...
        + sill3 * exp(-sqrt(Xlag.^2 + (alpha*Ylag).^2)/rho3));
elseif strcmp(model_2.cov1,'SumMetric')&&strcmp(model_2.cov2,'sph')
    sill3 = table2array(model_2.GeoVal(3,2));
    rho3 = table2array(model_2.GeoVal(3,3));
    alpha = table2array(model_2.GeoVal(3,4));
    Zlag = sqrt(Xlag.^2 + (alpha*Ylag).^2);
    fitVariog = nugget + sill1*(1.5*Xlag/rho1 - 0.5*(Xlag/rho1).^3)...
        + sill2*(1.5*Ylag/rho2 - 0.5*(Ylag/rho2).^3)...
        + sill3*(1.5*Zlag/rho2 - 0.5*(Zlag/rho2).^3);
    fitVariog(Xlag>rho1 & Ylag<=rho2 & Zlag<=rho3) = nugget + sill1...
        + sill2*(1.5*Ylag(Xlag>rho1 & Ylag<=rho2 & Zlag<=rho3)/rho2 - 0.5*(Ylag(Xlag>rho1 & Ylag<=rho2 & Zlag<=rho3)/rho2).^3)...
        + sill3*(1.5*Zlag(Xlag>rho1 & Ylag<=rho2 & Zlag<=rho3)/rho3 - 0.5*(Zlag(Xlag>rho1 & Ylag<=rho2 & Zlag<=rho3)/rho3).^3);
    fitVariog(Xlag>rho1 & Ylag>rho2 & Zlag<=rho3) = nugget + sill1 + sill2...
        + sill3*(1.5*Zlag(Xlag>rho1 & Ylag>rho2 & Zlag<=rho3)/rho3 - 0.5*(Zlag(Xlag>rho1 & Ylag>rho2 & Zlag<=rho3)/rho3).^3);
    fitVariog(Xlag>rho1 & Ylag>rho2 & Zlag>rho3) = nugget + sill1 + sill2 + sill3;
    fitVariog(Xlag<=rho1 & Ylag>rho2 & Zlag<=rho3) = nugget...
        + sill1*(1.5*Xlag(Xlag<=rho1 & Ylag>rho2 & Zlag<=rho3)/rho1 - 0.5*(Xlag(Xlag<=rho1 & Ylag>rho2 & Zlag<=rho3)/rho1).^3)...
        + sill2 + sill3*(1.5*Zlag(Xlag<=rho1 & Ylag>rho2 & Zlag<=rho3)/rho2 - 0.5*(Zlag(Xlag<=rho1 & Ylag>rho2 & Zlag<=rho3)/rho2).^3);
    fitVariog(Xlag<=rho1 & Ylag>rho2 & Zlag>rho3) = nugget...
        + sill1*(1.5*Xlag(Xlag<=rho1 & Ylag>rho2 & Zlag>rho3)/rho1 - 0.5*(Xlag(Xlag<=rho1 & Ylag>rho2 & Zlag>rho3)/rho1).^3)...
        + sill2 + sill3;
    fitVariog(Xlag<=rho1 & Ylag<=rho2 & Zlag>rho3) = nugget...
        + sill1*(1.5*Xlag(Xlag<=rho1 & Ylag<=rho2 & Zlag>rho3)/rho1 - 0.5*(Xlag(Xlag<=rho1 & Ylag<=rho2 & Zlag>rho3)/rho1).^3)...
        + sill2*(1.5*Ylag(Xlag<=rho1 & Ylag<=rho2 & Zlag>rho3)/rho2 - 0.5*(Ylag(Xlag<=rho1 & Ylag<=rho2 & Zlag>rho3)/rho2).^3)...
        + sill3;
% elseif strcmp(model_2.cov1,'ProdSum')&&strcmp(model_2.cov2,'matern')
%     k = table2array(model_2.GeoVal(1,5));
%     nu1 = table2array(model_2.GeoVal(1,4));
%     nu2 = table2array(model_2.GeoVal(2,4));
%     matern1 = 1/((2^(nu1-1))*gamma(nu1)) * ((2*sqrt(nu1)*Xlag)/rho1).^nu1 .* besselk(nu1,(2*sqrt(nu1)*Xlag)/rho1);
%     matern2 = 1/((2^(nu2-1))*gamma(nu2)) * ((2*sqrt(nu2)*Ylag)/rho2).^nu2 .* besselk(nu2,(2*sqrt(nu2)*Ylag)/rho2);
%     fitVariog = nugget+k*sill1*sill2+sill1+sill2...
%         -(sill1*matern1+sill2*matern2+k*(sill1*matern1).*(sill2*matern2));
% elseif strcmp(model_2.cov1,'ProdSum')&&strcmp(model_2.cov2,'exp')
%     k = table2array(model_2.GeoVal(1,4));
%     fitVariog = nugget + (k*sill2+1)*sill1*(1-exp(-Xlag/rho1))...
%     + (k*sill1+1)*sill2*(1-exp(-Ylag/rho2))...
%     - k*(sill1*(1-exp(-Xlag/rho1))).*(sill2*(1-exp(-Ylag/rho2)));
% elseif strcmp(model_2.cov1,'ProdSum')&&strcmp(model_2.cov2,'sph')
%     k = table2array(model_2.GeoVal(1,4)); % Please comment out when using expProdSum model
%     fitVariog = nugget + sill1*(1.5*Xlag/rho1 - 0.5*(Xlag/rho1).^3)...
%         + sill2*(1.5*Ylag/rho2 - 0.5*(Ylag/rho2).^3)...
%         + k*(sill1*(1.5*Xlag/rho1 - 0.5*(Xlag/rho1).^3)).*(sill2*(1.5*Ylag/rho2 - 0.5*(Ylag/rho2).^3));
%     fitVariog(Xlag>rho1 & Ylag<=rho2) = nugget + sill1...
%         + sill2*(1.5*Ylag(Xlag>rho1 & Ylag<=rho2)/rho2 - 0.5*(Ylag(Xlag>rho1 & Ylag<=rho2)/rho2).^3)...
%         + k*sill1.*(sill2*(1.5*Ylag(Xlag>rho1 & Ylag<=rho2)/rho2 - 0.5*(Ylag(Xlag>rho1 & Ylag<=rho2)/rho2).^3));
%     fitVariog(Xlag<=rho1 & Ylag>rho2) = nugget + sill1*(1.5*Xlag(Xlag<=rho1 & Ylag>rho2)/rho1 - 0.5*(Xlag(Xlag<=rho1 & Ylag>rho2)/rho1).^3)...
%         + sill2...
%         + k*(sill1*(1.5*Xlag(Xlag<=rho1 & Ylag>rho2)/rho1 - 0.5*(Xlag(Xlag<=rho1 & Ylag>rho2)/rho1).^3)).*(sill2);
%     fitVariog(Xlag>rho1 & Ylag>rho2) = nugget + sill1 + sill2 + k*sill1*sill2;
end
    
figure;
mesh(xlag,ylag,fitVariog,'FaceLighting','gouraud','LineWidth',0.5)
title('Fitted variogram')
xlabel('{\it x} lag (m)')
ylabel('{\it y} lag (m)')
zlabel('Semi-variance (t ha^{-1})^{2}')
