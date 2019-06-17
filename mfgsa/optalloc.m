% determines number of model evaluations and weights for multifidelity
% variance and sensitivity index estimation

% INPUTS
% p         computational budget (divide by (d+2) for effective budget for
%           Sobol sensitivity analysis)
% w         k-by-1 vector of computational costs for each model
% stats     struct containing model statistics
% varopt    if true, solve constrained minimization to yield m, alpha that
%           minimize variance of variance estimator. if false, return
%           m, alpha that minimize variance of mean estimator. false by
%           default.

% OUTPUTS
% m         k-by-1 vector of number of times to evaluate each model
% alpha     k-by-1 vector of multifidelity estimator weights

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 17 June 2019

function [m,alpha] = optalloc(p,w,stats,varopt)

if nargin == 3
    varopt = false;    
end

k = length(w);

% mean-optimal allocation (used as initial guess for variance optimization)
% see   B. Peherstorfer & K. Willcox, Optimal Model Management for
%       Multifidelity Monte Carlo Estimation, SISC 2016
temp = (stats.rho.^2 - [stats.rho(2:end); 0].^2);
r = sqrt(w(1)*temp./(w*(1-stats.rho(2)^2)));
m1 = p/(w'*r);
m  = floor(m1*r);
alpha = stats.rho*stats.sigma(1)./stats.sigma;

delta = stats.delta;
tau = stats.tau;
q = stats.q;
sigma = stats.sigma;
rho = stats.rho;

if varopt % variance-optimal allocation
    obj = @(m,a) 1/m(1)*(delta(1) - (m(1)-3)/(m(1)-1)*sigma(1)^4) + ...
            sum(a(2:end).^2.*(1./m(1:end-1).*(delta(2:end) - (m(1:end-1)-3)./(m(1:end-1)-1).*sigma(2:end).^4) ...
                - 1./m(2:end).*(delta(2:end) - (m(2:end)-3)./(m(2:end)-1).*sigma(2:end).^4))) ...
        + 2*sum(a(2:end).*...
            (1./m(2:end).*(q(2:end).*tau(1).*tau(2:end) + 2./(m(2:end)-1).*rho(2:end).^2.*sigma(1)^2.*sigma(2:end).^2) ...
           - 1./m(1:end-1).*(q(2:end).*tau(1).*tau(2:end) + 2./(m(1:end-1)-1).*rho(2:end).^2.*sigma(1)^2.*sigma(2:end).^2)));
    
    f = @(x) obj(x(1:k),x(k+1:end));
    
    % equality constraints to use up budget and to set alpha_1 = 1
    Aeq = [w', zeros(1,k); ...
           zeros(1,k), 1, zeros(1,k-1)];
    Beq = [p; 1];
    
    % inequality constraints so that m_{i+1} > m_i
    A = [[diag(ones(k-1,1)), zeros(k-1,1)] - [zeros(k-1,1), diag(ones(k-1,1))],zeros(k-1,k)];
    B = zeros(k-1,1);
    
    % inequality constraint so that m_1 > 3 (for a somewhat sensible
    % variance calculation)
    A = [A; -1, zeros(1,2*k-1)];
    B = [B; -3];
    x = fmincon(f,[m;alpha],A,B,Aeq,Beq,zeros(2*k,1),[]);
    m = floor(x(1:k));
    alpha = x(k+1:end);
end
