function res = lineofbestfit_fract_exp(xdata, ydata, num, dispmin, dispmax)
% dispmin/dispmax - displacement of the minimum/ maximum value compared to 
% the data set. If you only want data to begin and end at the min/max input
% points, these should be =0
% num - the number of data points you want calculated for your line of best
% fit. This should match the number of data points of your x-axis for
% plotting

model = @(b, x) b(1) * (1 ./ x).^b(2);
xdata = xdata(~isnan(xdata)) ; 
ydata = ydata(~isnan(ydata)) ; 

% Initial guess for k and n
initial_guess = [1, 1]; % Starting values for k and n

% Fit the model to the data
options = optimset('Display', 'off'); % Suppress output
[params, resnorm] = lsqcurvefit(model, initial_guess, xdata, ydata, [], [], options);

% Display the fitted parameters
k = params(1);
n = params(2);
disp(['k = ', num2str(k), ', n = ', num2str(n)]);

xaxis = linspace(min(xdata)-dispmin, max(xdata)+dispmax, num) ; 

% table = ones(length(xaxis)) ; 
res = k * (1./ xaxis).^n ;