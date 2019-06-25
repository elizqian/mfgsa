function Z = generate_inputs(N)
% generates random inputs for Ishigami function example
%
% INPUT
% N     number of inputs to generate
%
% OUTPUT
% Z     N-by-3 matrix of random Ishigami function inputs
%
% AUTHOR
% Elizabeth Qian (elizqian@mit.edu) 14 June 2019

Z = rand(N,3)*2*pi-pi;