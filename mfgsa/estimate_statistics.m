% estimate statistics of high- and low-fidelity models using N samples

% INPUTS
% fcns      cell array of anonymous functions corresponding to high- and
%           low-fidelity models (hi-fid must be first)
% N         number of samples to use to estimate statistics

% OUTPUT
% stats     struct with fields {mu, sigma, delta, tau, q, rho}

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 14 June 2019

function stats = estimate_statistics(fcns,N)

Z = generate_inputs(N);

k = length(fcns);

f_vals = zeros(N,k);

% loop through models
for i = 1:k
    f_vals(:,i) = fcns{i}(Z);   % assumes vectorized functions; need to loop thru inputs if not vectorized
end

stats.mu      = mean(f_vals)';
stats.sigma   = std(f_vals)';

g_vals        = (f_vals - stats.mu).^2;

a             = N^2/((N-1)*(N^2-3*N+3));
b             = 3*(2*N-3)/(N^2-3*N+3);
stats.delta   = (a*sum(g_vals.^2) - b*stats.sigma.^4)';

stats.tau     = std(g_vals)';

stats.q       = [1; zeros(k-1,1)];
stats.rho     = [1; zeros(k-1,1)];
for i = 2:k
    stats.rho(i) = sum((f_vals(:,1)-stats.mu(1)).*(f_vals(:,i)-stats.mu(i)))/((N-1)*stats.sigma(1)*stats.sigma(i));
    stats.q(i)   = (g_vals(:,1)-stats.sigma(1)^2)'*(g_vals(:,i)-stats.sigma(i)^2)/((N-1)*stats.tau(1)*stats.tau(i));
end





