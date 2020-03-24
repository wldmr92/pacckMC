%% Calculate critical points based on hopping channels
%  input parameters: lattice -> lattice object
%                    lut -> lookup table object for rate catalog
%                    xi -> minimal critical process propability for defining superbasins
%                    particle -> electron/hole
%  output parameters: lattice -> updated lattice object
function lattice = calc_critical_points(lattice,lut,xi,particle)

c = pp_constants;
% transitions and number of transitions
trans = {lut.trans}; n_trans = numel(trans);
% copy of lattice object
lattice_c = lattice;

for j = 1:1:c.x_size
    tic
    for i = 1:1:c.y_size
        pos = [i j];
        %% starting node i
        [~,idx_c] = max(lattice(pos(1),pos(2)).src_f(:));
        %% adjacent node j
        pos_t  = calc_transition(pos,trans{idx_c});
        % search only in one material (== exclude grain boundaries)
     %   if(strcmp(lattice(pos(1),pos(2)).struc,lattice(pos_t(1),pos_t(2)).struc) == 0)
     %       continue;
     %   else
            % index of auxiliary node vector
            n = 2;
            % store current hopping channel in auxiliary node vector
            lv(1) = lattice_c(pos(1),pos(2));
            lv(2) = lattice_c(pos_t(1),pos_t(2));
            % mark node i as local critical points
            lv(1).lcp = 1;
            % mark escape branches of node i
            lv(1).lcb(:) = 1;
            % mark critical branch of node i
            lv(1).lcb(idx_c) = 2;
            % mark node j as local critical points
            lv(2).lcp = 1;
            % mark escape branches of node j
            lv(2).lcb(:) = 1;
            % mark critical branch of node j
            idx_c_t = calc_reverse_idx(idx_c);
            lv(2).lcb(idx_c_t) = 2;
            % calculate critical process probability
            [~,p_c,f] = calc_probabilities(lv,n_trans,'l',particle);
            if p_c > xi
                % expand hopping channel by maximal escape flux
                [lattice] = add_maximal_escape_flux(lv,lattice,lattice_c,xi,p_c,n,trans,f,particle);
            else
                continue;
            end
      %  end
    end
    toc
end

end


 

             
                    
                    



