%% Calculate the scaling factor of the charge carrier within the superbasin
%  input parameter: alpha_sb -> maximum scaling factor of superbasin
%                   msbi -> number of scalings within this superbasin
%                   kappa -> scaling steps
%  output parameter: alpha -> scaling factor for charge carrier

function alpha = update_cc_alpha(alpha_sb,msbi,kappa)

    alpha = min(alpha_sb, alpha_sb^(msbi/kappa));

end