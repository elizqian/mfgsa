% computes multifidelity estimate of mean, variance, and Sobol' main and
% total effect sensitivity indices

% INPUTS
% fcns      k-by-1 cell array of functions corresponding to different models
% d         dimension of uncertaint input
% w         k-by-1 vector of computational costs for functions in fcns
% stats     struct containing statistics of models in fcns
% p         total computational budget
% vec       k-by-1 boolean vector that indicates whether the models in fcns
%           are vectorized or not. False is default.

% OUTPUTS
% mu        multifidelity mean estimate of high-fidelity model in fcns{1}
% sigsq     multifid. variance estimate of high-fidelity model in fcns{1}
% sm        d-by-1 vector of Sobol' main effect sensitivity indices
% st        d-by-1 vector of Sobol' total effect sensitivity indices

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 17 June 2019

function [mu,sigsq,sm,st] = mfsobol(fcns,d,w,stats,p,vec)

if nargin == 5
    vec = zeros(size(fcns));
end

% get optimal number of evaluations and weights using effective budget
[m,alpha] = optalloc(p/(d+2),w,stats);

% get two sets of independent inputs (generate_inputs.m for problem must be
% on MATLAB search path)
Z_A = generate_inputs(m(end));    
Z_B = generate_inputs(m(end));

% compute all evaluations of high-fidelity model
% this code assumes all models are vectorized - if models are not
% vectorized, need rewrite code to loop through inputs
if vec(1)
    yA = fcns{1}(Z_A(1:m(1),:));
    yB = fcns{1}(Z_B(1:m(1),:));
else
    yA = zeros(m(1),1); 
    yB = zeros(m(1),1);
    for j = 1:m(1)
        yA(j) = fcns{1}(Z_A(j,:));
        yB(j) = fcns{1}(Z_B(j,:));
    end
end

yC = zeros(m(1),d);
for i = 1:d
    Z_Ci = Z_B(1:m(1),:);
    Z_Ci(:,i) = Z_A(1:m(1),i);
    
    if vec(1)
        yC(:,i) = fcns{1}(Z_Ci);
    else
        for j = 1:m(1)
            yC(j,i) = fcns{1}(Z_Ci(j,:));
        end
    end
end

% initialize all statistics with their high-fidelity values
mu      = mean([yA; yB]);
sigsq   = var([yA; yB]);
[sm,st] = estimate_sobol(yA,yB,yC);

% loop through low-fidelity models
for j = 2:length(m)
    
    % get function evalutions - again, this code is for vectorized
    % functions
    if vec(j)
        yA = fcns{j}(Z_A(1:m(j),:));
        yB = fcns{j}(Z_B(1:m(j),:));
    else
        yA = zeros(m(j),1); 
        yB = zeros(m(j),1);
        for k = 1:m(j)
            yA(k) = fcns{j}(Z_A(k,:));
            yB(k) = fcns{j}(Z_B(k,:));
        end
    end

    yC = zeros(m(j),d);
    for i = 1:d
        Z_Ci = Z_B(1:m(j),:);
        Z_Ci(:,i) = Z_A(1:m(j),i);
        
        if vec(j)
            yC(:,i) = fcns{j}(Z_Ci);
        else
            for k = 1:m(j)
                yC(k,i) = fcns{j}(Z_Ci(k,:));
            end
        end
    end
    
    % add low-fi correction to existing estimate
    mu    = mu + alpha(j)*(mean([yA; yB]) - mean([yA(1:m(j-1)); yB(1:m(j-1))]));
    sigsq = sigsq + alpha(j)*(var([yA; yB]) - var([yA(1:m(j-1)); yB(1:m(j-1))]));
    
    [sm1,st1] = estimate_sobol(yA,yB,yC);
    [sm2,st2] = estimate_sobol(yA(1:m(j-1)),yB(1:m(j-1)), yC(1:m(j-1),:));
    
    sm    = sm + alpha(j)*(sm1-sm2);
    st    = st + alpha(j)*(st1-st2);
end