%% Check nearest neighbors for critical points and distribute superbasin indices
%  input parameters:  lattice -> lattice object
%                     pos -> current position of charge carrier
%                     trans -> possible hopping transitions
%  output parameters: lattice -> updated lattice object
function lattice = check_nearest_neighbors(lattice,pos,trans)

c = pp_constants;
% get number of local superbasins
idx_l = lattice(pos(1),pos(2)).lsbi;
% get number of global superbasins
idx_g = lattice(pos(1),pos(2)).gsbi;

% iterate over all possible transitions
for k = 1:4      
    % get new position
    pos_t = calc_transition(pos,trans{k});
    % check local superbasins
    if idx_l > 0
        if( pos(1) == pos_t(1) && pos(2) == pos_t(2) )
            continue;
        elseif(lattice(pos(1),pos(2)).lcp == 1 && lattice(pos_t(1),pos_t(2)).lcp == 1)
            lattice(pos(1),pos(2)).lcb(k) = 2;
        elseif(lattice(pos(1),pos(2)).lcp == 1 && lattice(pos_t(1),pos_t(2)).lcp == 0)
            lattice(pos(1),pos(2)).lcb(k) = 1;
            kp = calc_reverse_idx(k);
            lattice(pos_t(1),pos_t(2)).lcb(kp) = 1;
        end
        % group neighboring local critical points as local superbasin
        if (lattice(pos(1),pos(2)).lcp == 1 && lattice(pos_t(1),pos_t(2)).lcp == 1 && lattice(pos_t(1),pos_t(2)).lsbi == 0)
            lattice(pos_t(1),pos_t(2)).lsbi = lattice(pos(1),pos(2)).lsbi;
            % recursive call of function
            lattice = check_nearest_neighbors(lattice,pos_t,trans);
        end
    end
    % check global superbasins
    if idx_g > 0
        if( pos(1) == pos_t(1) && pos(2) == pos_t(2) )
            continue;
        elseif(lattice(pos(1),pos(2)).gcp == 1 && lattice(pos_t(1),pos_t(2)).gcp == 1)
            lattice(pos(1),pos(2)).gcb(k) = 2;
        elseif(lattice(pos(1),pos(2)).gcp == 1 && lattice(pos_t(1),pos_t(2)).gcp == 0)
            lattice(pos(1),pos(2)).gcb(k) = 1;
            kp = calc_reverse_idx(k);
            lattice(pos_t(1),pos_t(2)).gcb(kp) = 1;
        end
        % group neighboring global critical points as global superbasin
        if (lattice(pos(1),pos(2)).gcp == 1 && lattice(pos_t(1),pos_t(2)).gcp == 1 && lattice(pos_t(1),pos_t(2)).gsbi == 0)
            lattice(pos_t(1),pos_t(2)).gsbi = lattice(pos(1),pos(2)).gsbi;
            % recursive call of function
            lattice = check_nearest_neighbors(lattice,pos_t,trans);
        end
        
    end
end
end