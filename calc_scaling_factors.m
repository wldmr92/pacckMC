%% Calculate scaling factors and number of critical sightings for superbasins
%  input parameters: lattice -> lattice object
%                    idx_g -> number of global superbasins
%                    idx_l -> number of local superbasins
%                    gamma -> minimal critical process propability for defining superbasins
%                    delta -> error measure 
%  output parameters: lattice -> updated lattice object
function lattice = calc_scaling_factors(lattice,idx_g,idx_l,xi,delta,kappa,particle)

c = pp_constants;
n_trans = 4;

% interate over all local superbasins
for i = 1:idx_l
    % get all static rates in superbasin
    aux = reshape([lattice(:).lsbi] == i,c.y_size,c.x_size); [row,col] = find(aux == 1);
    % calculate escape probability
    [~,p_c] = calc_probabilities([lattice(aux ~= 0)],n_trans,'l',particle);        
    if p_c > xi
        alpha = calc_alpha(delta,p_c);
        N_f = calc_N_f(delta,alpha);
        for j = 1:numel(row)
            lattice(row(j),col(j)).alpha_l = alpha;
            lattice(row(j),col(j)).N_fl    = ceil(N_f);
            lattice(row(j),col(j)).p_c_l   = p_c;
            lattice(row(j),col(j)).p_es_l  = 1 - p_c;           
        end
    end 
end

% interate over all global superbasins
for i = 1:idx_g
    % get all static rates in superbasin
    aux = reshape([lattice(:).gsbi] == i,c.y_size,c.x_size); [row,col] = find(aux == 1);
    % calculate escape probability
    [~,p_c] = calc_probabilities([lattice(aux ~= 0)],n_trans,'g',particle);         
    if p_c > xi
        alpha = calc_alpha(delta,p_c);
        N_f = calc_N_f(delta,alpha);
        for j = 1:numel(row)
            lattice(row(j),col(j)).alpha_g = alpha;
            lattice(row(j),col(j)).N_fg    = ceil(N_f);
            lattice(row(j),col(j)).p_c_g   = p_c;
            lattice(row(j),col(j)).p_es_g  = 1 - p_c;              
        end
    end 
end

end