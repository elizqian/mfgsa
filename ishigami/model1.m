% high-fidelity model for Ishigami function example
% f1 = sin(z1) + a*sin(z2)^2 + b*z3^4*sin(z1)

% INPUTS
% Z         N-by-3 matrix of uncertain parameters, distributed ~U[-pi,pi]
% a,b       (optional) Ishigami function parameters, default a = 5, b = 0.1

% OUTPUT
% f1        N-by-1 vector of model evaluations at the uncertain inputs in Z

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 14 June 2019

function f1 = model1(Z,a,b)

if nargin == 1
    a = 5;
    b = 0.1;
end

f1 = sin(Z(:,1)) +  a*sin(Z(:,2)).^2 + b*Z(:,3).^4.*sin(Z(:,1));