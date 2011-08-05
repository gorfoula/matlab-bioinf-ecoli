function [SAMPLES] = RandLevels(n,m) %% generate uniformly distributed levels 0:m-1
%% INPUT
%%  n: number of samples
%%  m: number of different levels
SAMPLES=ceil(m*rand(n,2))-1;  %% generate random integers 0:m-1 from uniform distribution
end