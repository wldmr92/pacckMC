%% Calculate probability of residing in state i of superbasin B
%  input parameters: F_i -> free energy of state i
%                    Z -> partition function of superbasin B
%  output parameters: p_i -> probability of residing in state i 
function p_i = calc_state_probability(F_i,Z,particle)

c = pp_constants;
% thermal energy
E_th = c.k_B.*c.T;
% state probability
if(strcmp(particle,'electron') == 1)
    p_i = exp(-F_i./E_th)./Z;    
elseif(strcmp(particle,'hole') == 1) 
    p_i = 1 - exp(-F_i./E_th)./Z;    
else
   warning('Unknown particle!!!'); return; 
end

end