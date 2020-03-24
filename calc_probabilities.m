%% Calculate escape and critical process probability
%  input parameters: lv -> nodes of possible superbasin
%                    n_trans -> number of possible hopping transitions
%                    s -> 'l' == local or 'g'== global superbasin
%  output parameters: p_es -> escape process probability
%                     p_c  -> critical process probability
function [p_es,p_c,f] = calc_probabilities(lv,n_trans,s,particle)

% get static rate constants
src = [lv.src_f];
if strcmp(s,'l') == 1
    % get local critical branches
    cb = [lv.lcb];
elseif strcmp(s,'g') == 1
    % get global critical branches
    cb = [lv.gcb];
else
    warning('Unknown superbasin type!');
end
% store free energies
F = [lv.pes];
% calculate partition function
Z = calc_partition_function(F);
% calculate state probabilities
Pi = repelem(calc_state_probability(F,Z,particle),n_trans);
% calculate probability fluxes
f = src.*Pi;
% calculate escape event probability flux
f_es = sum(f(cb == 1));
% calculate critical event probability flux
f_c  = sum(f(cb == 2));
% escape process probability 
p_es = f_es/(f_c + f_es);
% critical process probability
p_c = f_c/(f_c + f_es);
    
end