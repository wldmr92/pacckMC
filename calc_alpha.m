%% Calculate the scaling factor of a superbasin
%  input parameter: delta -> error measure
%                   p_c -> critical process probability
%  output parameter: alpha -> scaling factor

function alpha = calc_alpha(delta,p_c)

    alpha = ((delta./(1-delta)).*(p_c./(1-p_c)) + 1./(1-delta));

end