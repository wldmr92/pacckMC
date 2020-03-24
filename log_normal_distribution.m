%% Calculation of characteristic parameters of a log-normal distribution
%  input parameters: mu_L -> mean of transformed random variable Y = log(X)
%                    s_L  -> standard deviation of transformed random variable Y = log(X)
%  output parameters: mu -> mean of log-normal distribution of random variable X
%                     s  -> standard deviation of log-normal distributed of random variable X
function [mu,s] = log_normal_distribution(mu_L,s_L)

% mean of log-normal distribution
mu = exp(mu_L + ((s_L.^2)/2));
% standard deviation of log-normal distribution
s = (exp(s_L.^2) - 1).*exp(2.*mu_L + (s_L.^2));

end