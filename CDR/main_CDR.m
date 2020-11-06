% computes multifidelity sensitivity index estimates for
% convection-diffusion-reaction example, bootstrapping from 120000
% precomputed (yA,yB,yC) samples

% PAPER
% E. Qian, B. Peherstorfer, D. O'Malley, V. Vesselinov, and K. Willcox
% Multifidelity estimation of variance and sensitivity indices
% SIAM/ASA Journal on Uncertainty Quantification, 6(2):683-706, 2018.

% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 17 June 2019

%% SETUP
clear
addpath('../mfgsa')
samples = load('samples.mat');  % load pre-computed samples for bootstrapping

d = 5;                          % dimension of uncertain input

% function definitions that bootstrap from precomputed function outputs
fcns{1} = @(Z) deal(samples.yA(Z,1), samples.yB(Z,1), squeeze(samples.yC(Z,1,:)));
fcns{2} = @(Z) deal(samples.yA(Z,2), samples.yB(Z,2), squeeze(samples.yC(Z,2,:)));

w   = [1.94; 6.2e-3];       % assign model weights/costs
vec = [2 2];                % says that functions are bootstrapping

budget      = 1000*60;       % minutes times seconds

%% COMPUTE MULTIFIDELITY GLOBAL SENSITIVITY ANALYSIS REPLICATES
n_reps = 100;    % number of replicates 
estim  = 'Saltelli'; % which estimator to use -- 'Owen' or 'Saltelli'

% allocate storage
avg   = zeros(n_reps,2);    vr    = zeros(n_reps,2);
mc_sm = zeros(n_reps,d);    mc_st = zeros(n_reps,d);
mf_sm = zeros(n_reps,d);    mf_st = zeros(n_reps,d);

for n = 1:n_reps
    % estimate model statistics using small pilot sample
    stats = estimate_statistics(fcns,10,vec);
    
    % call mfsobol.m with just the high-fidelity model to get Monte
    % Carlo estimate
    [sm,st,mu,sigsq] = mfsobol(fcns(1),d,w(1),stats,budget,vec(1),estim);
    avg(n,1) = mu; 
    vr(n,1) = sigsq;
    mc_sm(n,:) = sm;
    mc_st(n,:) = st;
    
    % call mfsobol.m with full array of functions to get multifidelity
    % estimates
    [sm,st,mu,sigsq] = mfsobol(fcns,d,w,stats,budget,vec,estim);
    avg(n,2) = mu; 
    vr(n,2) = sigsq;
    mf_sm(n,:) = sm;
    mf_st(n,:) = st;
end

%% BOXPLOTS
warning('off','MATLAB:legend:IgnoringExtraEntries')
blue = [0       0.4470 0.7410];
red  = [0.8500  0.3250 0.0908];

% plot main effect sensitivity indices
figure(1); clf
plot(0,0,'Color',red); hold on; plot(0,0,'Color',blue);
h = boxplot([mc_sm(:,1), mc_sm(:,2), mc_sm(:,3), mc_sm(:,4), mc_sm(:,5), mf_sm(:,1), mf_sm(:,2), mf_sm(:,3), mf_sm(:,4), mf_sm(:,5)],...
    'Colors',[red; red; red; red; red; blue; blue; blue; blue; blue],'Whisker',10,...
    'labels',{'$s_m^1$','$s_m^2$','$s_m^3$','$s_m^4$','$s_m^5$','$s_m^1$','$s_m^2$','$s_m^3$','$s_m^4$','$s_m^5$'});
set(h,{'linew'},{2}); grid on
legend({'Monte Carlo','Multifidelity '},...
    'Location','NorthWest','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title(['CDR main sensitivity indices: $p = $',num2str(budget/60),' min'],'interpreter','latex')

% plot total effect sensitivity indices
figure(2); clf
plot(0,0,'Color',red); hold on; plot(0,0,'Color',blue);
h = boxplot([mc_st(:,1), mc_st(:,2), mc_st(:,3), mc_st(:,4), mc_st(:,5), mf_st(:,1), mf_st(:,2), mf_st(:,3), mf_st(:,4), mf_st(:,5)],...
    'Colors',[red; red; red; red; red; blue; blue; blue; blue; blue],'Whisker',10,...
    'labels',{'$s_t^1$','$s_t^2$','$s_t^3$','$s_t^4$','$s_t^5$','$s_t^1$','$s_t^2$','$s_t^3$','$s_t^4$','$s_t^5$'});
set(h,{'linew'},{2}); grid on
legend({'Monte Carlo','Multifidelity '},...
    'Location','Best','interpreter','latex'); legend boxoff
bp = gca;   bp.XAxis.TickLabelInterpreter = 'latex';
title(['CDR total sensitivity indices, $p = $',num2str(budget/60),' min'],'interpreter','latex')
