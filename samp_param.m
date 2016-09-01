%%Sample from a prior: provide minimum and maximum guesses in column vector
%%form, and a total number of particles, N, to sample.

%uniform prior
function [theta] = samp_param(min,max,N)
dim = length(min);
% if dim ~= length(max)
%    error('dimensions of minimum guess and maximum guess do not match')
% end
diff = max-min;
min_N = repmat(min,1,N);
diff_N = repmat(diff,1,N);

theta = min_N + diff_N.*rand(dim,N);