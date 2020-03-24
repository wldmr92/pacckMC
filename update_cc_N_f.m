%% Update number of critical sightings for charge carrier
%  input parameter: delta -> error measure
%                   alpha_sb -> maximum scaling factor of SB
%                   msbi -> number of scalings already performed
%  output parameter: N_f -> number of critical sightings

function N_f = update_cc_N_f(delta,alpha_sb,msbi,kappa)

    N_f = log(1./delta).*(1 - (1./delta) + sqrt(update_cc_alpha(alpha_sb,msbi,kappa)).*((1-delta)./delta));
    %N_fg = lattice(pos(1),pos(2)).N_fg;
end 