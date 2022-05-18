function [sm,st] = estimate_sobol(method,yA,yB,X)
% INPUTS
% method    'Owen' or 'Saltelli' or 'Gamboa'
% yA        N-by-1 first set of function evaluations
% yB        for 'Owen' or 'Saltelli' only: N-by-1 second set of function 
%           evaluations. For 'Gamboa' pass an empty array []
% X         **for 'Owen' or Saltelli': N-by-d function evaluationsthe i-th 
%           column of X was evaluated at the same inputs that led to yB 
%           except replacing the i-th column of the input with the i-th 
%           column of yA. 
%           **for 'Gamboa': N-by-d array of inputs that led to outputs yA
%
% OUTPUTS
% sm        d-by-1 vector of main effect sensitivities
% st        d-by-1 vector of total effect sensitivities (zero output for
%           method = 'Gamboa')
%
% AUTHORS
% Elizabeth Qian (www.elizabethqian.com)
% Giuseppe Cataldo
%
% LAST UPDATED
% 18 May 2022

[N,d] = size(X);

muA  = mean(yA);
varA = var(yA);

switch method
    case 'Owen'
        sm = 2*N/(2*N-1)*(yA'*X/N - (muA+mean(X)).^2/4 + (varA+var(X))/(4*N))/varA;
        st = 1/(2*N)*sum((yB-X).^2)/varA;
    case 'Saltelli'
        sm = (1/(N-1)*yA'*X - muA^2)/varA;
        st = 1 - (1/(N-1)*yB'*X - muA^2)/varA;
    case 'Gamboa'
        % Get X ordering
        [~,px] = sort(X);
        % Get pi_j
        [~,pi_j] = sort(px);
        % Get N_j with shift
        argpiinv = mod(pi_j,N)+1;
        % Index in x of thing with this rank
        N_j = zeros(N,d);
        for i = 1:d
            N_j(:,i) = px(argpiinv(:,i),i);
        end
        % Get output corresponding to new rank
        y_Nj = yA(N_j);
        sm = (1/N*yA'*y_Nj - muA^2)/varA;
        st = zeros(1,d);
end
