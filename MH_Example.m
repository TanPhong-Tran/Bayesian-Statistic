close all; clear all;

%% --- Create density function ----
DoubleExp = @(x) (1/2 * exp(-abs(x))); % same lambda function in python
Unif = @(y,x,delta) (unifpdf(y, x-delta/2, x+delta/2));

% --- Setting up MH ---
N = 100000; burn = 1000;
delta = 1;
x = zeros(1, N);
x(1) = unifrnd(-delta/2, delta/2); 

% --- Perform MH ---
for i = 2:N
    y = x(i-1) + unifrnd(-delta/2, delta/2);
    u = unifrnd(0,1);
    rho = DoubleExp(y)/DoubleExp(x(i-1)) * Unif(y,x(i-1),delta)/Unif(x(i-1),y,delta);
    x(i) = x(i-1) + (u<=rho)*(y-x(i-1));
end

% --- Extract post burn-in ---
x = x(burn+1:N);

% --- Plotting
xPlot = linspace(min(x), max(x), 1000);
histogram(x, 'Normalization', 'pdf');
hold on; 
plot(xPlot, DoubleExp(xPlot), 'LineWidth', 2);
hold off;

%% --- Norm-Cauchy ---
% Formulate X|theta ~ N(0,1); theta ~ Cauchy(0,1); 
% f(theta|X) ~appr 1/(1+theta^2) exp(-1/2 * (theta)- 1)^2

fCauchy = @(theta) (1./(1+theta.^2) .* exp(-1./2*(theta-1).^2));
fNorm = @(y,theta,sig) (exp(-1./(2.*sig^2)*(y-theta).^2));

% --- Setting up MH ---
N = 100000; burn = 1000;
sig = 1;
x = zeros(1, N);
x(1) = randn(1); 

% --- Perform MH ---
for i = 2:N
    y = x(i-1) + sig*randn(1);
    u = rand;
    rho = fCauchy(y)/fCauchy(x(i-1)) * fNorm(y,x(i-1),sig)/fNorm(x(i-1),y,sig);
    x(i) = x(i-1) + (u<=rho)*(y-x(i-1));
end

% --- Extract post burn-in ---
x = x(burn+1:N);
M = integral(fCauchy, -Inf, Inf); % calculate normalize term for integral of pdf

% --- Plotting ---
xPlot = linspace(min(x), max(x), 1000);
histogram(x, 'Normalization', 'pdf');
hold on; 
plot(xPlot, fCauchy(xPlot)/M, 'LineWidth', 2);
hold off;