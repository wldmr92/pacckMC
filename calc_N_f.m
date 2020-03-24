%% Calculate number of critical sightings
%  input parameter: delta -> error measure
%                   alpha -> scaling factor
%  output parameter: N_f -> number of critical sightings

function N_f = calc_N_f(delta,alpha)

    N_f = log(1./delta).*(1 - (1./delta) + alpha.*((1-delta)./delta));
    
end 