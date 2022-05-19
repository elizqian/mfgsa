% computes multifidelity sensitivity index and variance estimates for the
% Ishigami function and compares to the high-fidelity Monte Carlo estimate

% PAPERS
% E. Qian, B. Peherstorfer, D. O'Malley, V. Vesselinov, and K. Willcox
% Multifidelity estimation of variance and sensitivity indices
% SIAM/ASA Journal on Uncertainty Quantification, 6(2):683-706, 2018.
%
% G. Cataldo, E. Qian, and J. Auclair
% Multifidelity uncertainty quantification and model validation of
% large-scale multidisciplinary systems
% Journal of Astronomical Telescopes, Instruments and Systems, 2022
%
% AUTHORS
% Elizabeth Qian (elizqian@alum.mit.edu) 
% Giuseppe Cataldo
%
% LAST UPDATED
% 18 May 2022

%% SETUP
clear; close all
addpath('../mfgsa')

% define high-fidelity and low-fidelity models
fcns{1} = @(Z) model1(Z);   % high-fidelity
fcns{2} = @(Z) model2(Z);   % low-fidelity
fcns{3} = @(Z) model3(Z);   % lowest-fidelity

w = [1; 0.05; 0.001];       % assign model weights/costs
vec = [1 1 1];              % says each model is vectorized

d = 3;          % dimension of uncertain input

budget = 2000;   % define computational budget

% if analytical statistics are not available, estimate them. otherwise load
% true statistics from file
estimate = false;   
if estimate
    n_estimate = 100;
    stats = estimate_statistics(fcns,n_estimate);
else
    load('ishigami/truestats.mat');
end


%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 100;    % number of replicates 
avg   = zeros(n_reps,2);    vr    = zeros(n_reps,2);
mc_sm = zeros(n_reps,d);    mc_st = zeros(n_reps,d);
mf_sm = zeros(n_reps,d);    mf_st = zeros(n_reps,d);

method = 'Gamboa';

for n = 1:n_reps
    
    % call mfsobol.m with just the high-fidelity model to get Monte
    % Carlo estimate
    [sm,st,mu,sigsq] = mfsobol(fcns(1),d,w(1),stats,budget,vec(1),method);
    avg(n,1) = mu; 
    vr(n,1) = sigsq;
    mc_sm(n,:) = sm;
    mc_st(n,:) = st;
    
    % call mfsobol.m with full array of functions to get multifidelity
    % estimates
    [sm,st,mu,sigsq] = mfsobol(fcns,d,w,stats,budget,vec,method);
    avg(n,2) = mu; 
    vr(n,2) = sigsq;
    mf_sm(n,:) = sm;
    mf_st(n,:) = st;
end

%% PLOT ESTIMATOR SPREAD
warning('off','MATLAB:legend:IgnoringExtraEntries')
blue = [0       0.4470 0.7410];
red  = [0.8500  0.3250 0.0908];

% plot main effect sensitivity indices
figure(1); clf
h = boxplot([mc_sm(:,1), mf_sm(:,1), mc_sm(:,2), mf_sm(:,2), mc_sm(:,3), mf_sm(:,3)],...
    'Colors',[blue; red; blue; red; blue; red],'Whisker',10,...
    'labels',{'MC $s_m^1$','MF $s_m^1$','MC $s_m^2$','MF $s_m^2$','MC $s_m^3$','MF $s_m^3$'});
set(h,{'linew'},{2}); grid on
legend(flipud(findall(gca,'Tag','Box')), {'High-fidelity estimate','Multifidelity estimate'},...
    'Location','SouthWest','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title([method,' main sensitivity estimates for Ishigami function'],'interpreter','latex')

if strcmp(method,'Owen') || strcmp(method,'Saltelli')
    % plot total effect sensitivity indices
    figure(2); clf
    h = boxplot([mc_st(:,1), mf_st(:,1), mc_st(:,2), mf_st(:,2), mc_st(:,3), mf_st(:,3)],...
        'Colors',[blue; red; blue; red; blue; red],'Whisker',10,...
        'labels',{'MC $s_t^1$','MF $s_t^1$','MC $s_t^2$','MF $s_t^2$','MC $s_t^3$','MF $s_t^3$'});
    set(h,{'linew'},{2}); grid on
    hLegend = legend(flipud(findall(gca,'Tag','Box')), {'High-fidelity estimate','Multifidelity estimate'},...
        'Location','NorthEast','interpreter','latex'); legend boxoff
    bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
    title([method,' total sensitivity estimates for Ishigami function'],'interpreter','latex')
end

% plot variance estimates
figure(3); clf
histogram(vr(:,1),12,'facecolor',blue); hold on
histogram(vr(:,2),12,'facecolor',red,'facealpha',1)
legend({'High-fidelity','Multifidelity'},'interpreter','latex'); legend boxoff
title('Variance estimates for Ishigami function','interpreter','latex')