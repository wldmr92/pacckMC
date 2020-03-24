%% Calculation of static rate constants
%  input parameters: lattice -> lattice object
%                    lut -> lookup table object for rate catalog
%                    E_a_ref_g -> global superbasin reference energy barrier
%                    scaling -> string for on/off switching of scaling algorithm
%  output parameters: lattice -> updated lattice object
%                     E_a -> equivalent activation energies
%                     GAMMA -> wave function overlap parameters
function [lattice,E_a,GAMMA] = calc_static_rates(lattice,lut,E_a_ref_g,scaling)

c = pp_constants;
% number of transitions
n_trans = numel(lut);
% equivalent activation energies 
E_a = zeros(n_trans*c.x_size*c.y_size,1);
% wave function overlap parameters 
GAMMA = zeros(n_trans*c.x_size*c.y_size-2.*c.x_size,1);

idx = 1; idx_G = 1;
for j = 1:c.x_size
    for i = 1:c.y_size
        % current position
        pos = [i j];
       
        for k = 1:n_trans % for each transition
            % current transition position
            pos_t = calc_transition(pos,lut(k).trans);
            if ( pos(1) == pos_t(1) && pos(2) == pos_t(2) )
                lattice(i,j).src_f(k) = 0; % forward rate
                lattice(i,j).src_b(k) = 0; % backward rate
            else
                % energy difference between adjacent sites i and j
                Delta_E_ij = lattice(pos_t(1),pos_t(2)).pes - lattice(pos(1),pos(2)).pes;
                % calculate prefactor between adjacent sites i and j
                a_ij = sqrt( lattice(pos(1),pos(2)).a*lattice(pos_t(1),pos_t(2)).a );
                % calculate coupling between adjacent sites i and j
                Gamma_ij = lattice(pos(1),pos(2)).Gamma + lattice(pos_t(1),pos_t(2)).Gamma;
                GAMMA(idx_G) = Gamma_ij; idx_G = idx_G + 1;
                % calculate hopping distance between adjacent sites i and j
                Delta_R_ij = calc_distance(pos,pos_t);
                % calculate and store static rate constant
                % forward transition
                k_ij = calc_hopping_rates(Delta_E_ij,'hole',a_ij,Gamma_ij,Delta_R_ij);
                lattice(i,j).src_f(k) = k_ij;
                % backward transition
                k_ji = calc_hopping_rates(-Delta_E_ij,'hole',a_ij,Gamma_ij,Delta_R_ij);
                lattice(i,j).src_b(k) = k_ji;
                if(strcmp(scaling,'scaling_on') == 1)
                    % calculate effective activation energy of Arrhenius equation
                    E_a(idx) = calc_equivalent_activation_energy(k_ij);
                    % mark global critical points and branches
                    if(E_a(idx) > E_a_ref_g)
                        lattice(pos(1),pos(2)).gcp = 1;
                        lattice(pos(1),pos(2)).gcb(1:n_trans) = 1; lattice(pos(1),pos(2)).gcb(k) = 2;
                    end
                    idx = idx + 1;
                elseif(strcmp(scaling,'scaling_off') == 1)
                    continue;
                else
                    warning('Unknown input parameter!!!')
                end
            end
        end
    end
end



end