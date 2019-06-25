function [sm,st] = estimate_sobol(yA,yB,yC,method)
% computes main and total effect sensitivities from matrices of function
% evaluations yA,yB,yCi
%
% INPUTS
% yA        M-by-1 first set of function evaluations
% yB        M-by-1 second set of function evaluations
% yC        M-by-d function evaluations; the i-th column of yC was evaluated
%           at the same inputs that led to yB except replacing the i-th
%           column of the input with the i-th column of yA
% method    'Owen' or 'Saltelli'; uses Owen estimators by default.
%
% OUTPUTS
% sm        d-by-1 vector of main effect sensitivities
% st        d-by-1 vector of total effect sensitivities
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 17 June 2019

if nargin == 3
    method = 'Owen';
end

N = length(yA);

muA = mean(yA);
muB = mean(yB);
varA = var(yA);
varB = var(yB);

switch method
    case 'Owen'
        sm = 2*N/(2*N-1)*(yA'*yC/N - (muA+muB)^2/4 + (varA+varB)/(4*N))/varA;
        st = 1/(2*N)*sum((yB-yC).^2)/varA;
    case 'Saltelli'
        sm = (1/(N-1)*yA'*yC - muA^2)/varA;
        st = 1 - (1/(N-1)*yB'*yC - muA^2)/varA;
end