clear all; close all;

%% ----- Settings -----
N = 100000; 
burn = 5000;
beta = 2;
x = [0.2, 0.1, 0.25];
alphas =  zeros(1,N);
etas = zeros(1,N);
alphas(1) = exprnd(1);
etas(1) = exprnd(1);

%% ----- Create density functions
post = @(alpha, eta) (alpha*eta)^3 * ...
    exp(-eta*(sum(x.^alpha) + beta)-alpha)*eta^(beta-1)*(prod(x.^(alpha-1)));

q = @(alpha_h, eta_h, alpha, eta) (1/(alpha*eta) * ...
    exp(-(alpha_h/alpha) - (eta_h/eta)));

%% --- Perform MH ---
for i=2:N
    alpha_t = exprnd(alphas(i-1));
    eta_t = exprnd(etas(i-1));
    rho = (post(alpha_t, eta_t) / post(alphas(i-1), etas(i-1))) * ...
        (q(alphas(i-1), etas(i-1), alpha_t, eta_t) /...
        q(alpha_t, eta_t, alphas(i-1), etas(i-1)));
    u = rand;
    alphas(i) = alphas(i-1) + (u<=rho)*(alpha_t - alphas(i-1));
    etas(i) = etas(i-1) + (u<=rho)*(eta_t - etas(i-1));
end

%% Extract post burn-inburn-in
alphas = alphas(burn+1:N);
etas = etas(burn+1:N);

%% Plot histogram
subplot(1,2,1);
histogram(alphas,'Normalization','pdf');

subplot(1,2,2);
histogram(etas,'Normalization','pdf');