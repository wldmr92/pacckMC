%% Add the maximal escape rate to a given set of nodes
%  input parameter: lv -> auxiliary node vector for storing a possible superbasin
%                   lattice   -> current lattice object
%                   lattice_c -> copy of the empty lattice object
%                   gamma -> minimal critical process propability for defining superbasins
%                   p_c -> critical process probability of hopping channel
%                   trans -> possible hopping transitions
%  output parameter: lattice -> updated lattice objected
function [lattice] = add_maximal_escape_flux(lv,lattice,lattice_c,xi,p_c,n,trans,f,particle)

n_trans = numel(trans);

% add maximal escape probability flux until p_c is lower than xi or p_c is decreasing
while(1)
    % find maximal escape flux
    f([lv(1:n).lcb] == 2) = 0;
    f_max = max(f([lv(1:n).lcb] == 1)); f_max = f_max(1);
    idx_max  = find(f == f_max); idx_max = idx_max(1);
    if mod(idx_max,n_trans) ~= 0
        n_max = floor(idx_max/n_trans) + 1;
        idx_max = mod(idx_max,n_trans);
    else
        n_max = idx_max/n_trans;
        idx_max = n_trans;
    end
    % get node position corresponding to maximal escape flux (pos_max)
    pos_max = lv(n_max).pos;
    % calculate position of new node added to potential superbasin
    pos_t = calc_transition(pos_max,trans{idx_max});
    % check if new node is already element of lv()
    if sum([lv(:).id] == lattice(pos_t(1),pos_t(2)).id) == 0
        % increment node index
        n = n + 1;
        % add new node to lv()
        lv(n) = lattice_c(pos_t(1),pos_t(2));
        % mark new node as local critical point
        if lv(n).lcp == 0
            lv(n).lcp = 1;
            % mark escape branches of new node
            lv(n).lcb(:) = 1;
        end
    end
    % calculate reverse index
    idx_c_t = calc_reverse_idx(idx_max);
    % mark critical branch of node i
    lv(n).lcb(idx_c_t) = 2;
    % mark critical branch of node j
    lv(n_max).lcb(idx_max) = 2;
    % store old escape probability
    p_c_p = p_c;
    % calculate critical process probability
    [~,p_c,f] = calc_probabilities(lv,n_trans,'l',particle);
    % calculate absolute change in critical process probability
    Delta_p_c = p_c - p_c_p;
    % check if critical process probability is higher than treshold value xi AND increasing
    if(Delta_p_c > 0)       
        continue;
    else
        % update lattice object
        if p_c_p > xi
            for k = 1:n-1
                pos = lv(k).pos;
                lattice(pos(1),pos(2)) = lv(k);
            end
            break;
        end
    end
end

end