classdef particles 
    
    properties 
        
        type   = 'hole';
        id     = 0;
        pos    = [0, 0];
        alpha_l  = 1;
        alpha_g  = 1;
        N_fl  = 0;
        N_fg  = 0;
        scale_g = 0;
        scale_l = 0;
        lsb = 0;
        gsb = 0;
        ctr_as_lsb = 0;
        ctr_as_gsb = 0;
        lsbi  = [0, 0, 0];
        gsbi  = [0, 0, 0];
        
    end
    
    methods
        % constructor method
        function particle = particles(type,id,N_fl,N_fg)
            if nargin ~= 0
                m = size(id,1);
                n = size(id,2);
                particle(m,n) = particle;
                for i = 1:m
                    for j = 1:n
                        particle(i,j).type  = type{i,j};
                        particle(i,j).id    = id(i,j);
                        particle(i,j).N_fl = N_fl;
                        particle(i,j).N_fg = N_fg;                       
                    end
                end
            end
        end        
    end
       
end