%% Group critical points to superbasins
%  input parameters: lut -> lookup table object for rate catalog
%                    lattice -> lattice object 
%  ouput parameters: lattice -> updated lattice object
%                    idx_g -> number of global superbasins
%                    idx_l -> number of local superbasins
function [lattice,idx_g,idx_l] = group_superbasins(lut,lattice)

c = pp_constants;

idx_g = 0; idx_l = 0;
for j = 1:c.x_size
    for i = 1:c.y_size
        % current position
        pos = [i j];  
        % local superbasins
        if (lattice(i,j).lsbi == 0 && lattice(i,j).lcp == 1)
            idx_l = idx_l + 1;                     
            lattice(i,j).lsbi = idx_l;
            lattice = check_nearest_neighbors(lattice,pos,{lut.trans});      
        end
        % global superbasins
        if (lattice(i,j).gsbi == 0 && lattice(i,j).gcp == 1)
            idx_g = idx_g + 1;                     
            lattice(i,j).gsbi = idx_g;
            lattice = check_nearest_neighbors(lattice,pos,{lut.trans});   
            % check if global superbasin is bigger than one node
            if(sum(sum(reshape([lattice(:).gsbi] == idx_g,c.y_size,c.x_size))) == 1)
                lattice(i,j).gsbi = 0;
                idx_g = idx_g - 1; 
            end            
        else
            continue;
        end                  
    end    
end

end