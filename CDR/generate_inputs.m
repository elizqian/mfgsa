function Z = generate_inputs(N)
% generates random indices for bootstrapping from pre-computed CDR samples
%
% INPUT
% N     number of inputs to generate
%
% OUTPUT
% ind   N-by-1 vector of random indices
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 5 November 2020

Z = randi(120000,N,1); % 120000 is the number of pre-computed CDR samples