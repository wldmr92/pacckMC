%% Calculate hopping rates based on Miller-Abrahams equation
%  input parameters: Delta_E_ij -> energy difference between nodes i and j
%                    particle   -> particle type ('electron' or 'hole')
%                    a_ij       -> prefactor between node i and j (a_ij = sqrt[a_i*a_j])
%                    Gamma_ij   -> positional coupling between node i and j (Gamma_ij = Gamma_i + Gamma_j)
%                    Delta_R_ij -> hopping distance between node i and j 
%  output parameters: k_ij -> hopping rate between node i and j
function k_ij = calc_hopping_rates(Delta_E_ij,particle,a_ij,Gamma_ij,Delta_R_ij)

c = pp_constants;

if(strcmp(particle,'electron') == 1)           
        if(Delta_E_ij <= 0)
            % downwards hopping is favored
            k_ij = a_ij.*exp(-Gamma_ij.*(Delta_R_ij./c.r_L));
        else
            % upwards hopping is punished
            k_ij = a_ij.*exp(-Gamma_ij.*(Delta_R_ij./c.r_L)).*exp((-Delta_E_ij)./(c.k_B.*c.T));
        end
elseif(strcmp(particle,'hole') == 1)
        if(Delta_E_ij >= 0)
            % upwards hopping is favored
            k_ij = a_ij.*exp(-Gamma_ij.*(Delta_R_ij./c.r_L));
        else
            % downwards hopping is punished
            k_ij = a_ij.*exp(-Gamma_ij.*(Delta_R_ij./c.r_L)).*exp((Delta_E_ij)./(c.k_B.*c.T));
        end
end

end