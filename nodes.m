classdef nodes
   
   properties 
       
       id        = 0;                   % unique identifier
       mat       = 'P3HT';				% material
       struc     = 'amorph';			% structure -> amorph/crystalline
	   pos       = [0, 0];				% position as 2D vector
       pes 	     = 0;					% static potential energy surface (PES)
       Gamma     = 0;                   % (-> off-diagonal disorder) (1)
       a         = 0;                   % prefactor per node (-> off-diagonal disorder)
	   src_f     = [0, 0, 0, 0];        % static rate constants (-> forward)
       E_a_f     = [0, 0, 0, 0];        % equivalent activation energies (-> forward)
       src_b     = [0, 0, 0, 0];        % static rate constants (-> backward)
       E_a_b     = [0, 0, 0, 0];        % equivalent activation energies (-> backward)
       
       % attributes for local superbasins
	   lcb 	     = [0, 0, 0, 0];        % flags for local critical branches  -> 0 == non-critical, 1 == escape-branch, 2 == critical branch, 3 == joining branch between two critical nodes but non-critical     
	   lcp		 = 0;                   % flag for local critical point  
	   lsbi		 = 0;					% local superbasin index  
       alpha_l   = 1;                   % local scaling factor
       N_fl      = 0;                   % number of local critical sightings
       t_lsbi    = 0;                   % local superbasin time
       p_c_l     = 0;
       p_es_l    = 0;
       
       % attributes for global superbasins (e.g. crystalline regions)
	   gcb 	     = [0, 0, 0, 0];        % flags for global critical branches -> 0 == non-critical, 1 == escape-branch, 2 == critical branch, 3 == joining branch between two critical nodes but non-critical        
	   gcp		 = 0;                   % flag for global critical point            
       gsbi      = 0;                   % global superbasin index   
       alpha_g   = 1;                   % global scaling factor 
       N_fg      = 0;                   % number of global critical sightings
       t_gsbi    = 0;                   % global superbasin time  
       p_c_g     = 0;
       p_es_g    = 0;       
   end
   
   methods
        % constructor method
        function node = nodes(id,material,structure,X,Y,PES,a,Gamma)
            if nargin ~= 0
                m = size(material,1);
                n = size(material,2);
                node(m,n) = node;
                for i = 1:m
                    for j = 1:n
                        node(i,j).id     = id(i,j);
                        node(i,j).mat    = material{i,j};
                        node(i,j).struc  = structure{i,j};
                        node(i,j).pos(1) = Y(i,j);
                        node(i,j).pos(2) = X(i,j);
                        node(i,j).pes    = PES(i,j);
                        node(i,j).a      = a(i,j);
                        node(i,j).Gamma  = Gamma(i,j);                        
                    end
                end
            end
        end 
        
        % display method 
        function h = display(nodes,X,Y)
            if nargin ~= 0
                m = size(nodes,1);
                n = size(nodes,2); 
                lcp = zeros(m,n);
                lsbi = zeros(m,n);
                lsbi_helper = zeros(m,n);
                gcp = zeros(m,n);
                gsbi = zeros(m,n);      
                p_c_l  = zeros(m,n);  
                p_es_l = zeros(m,n); 
                p_c_g = zeros(m,n); 
                p_es_g = zeros(m,n); 
                alpha = zeros(m,n); 
                for i = 1:m
                   for j = 1:n
                       pes(i,j)  = nodes(i,j).pes;
                       lcp(i,j)  = nodes(i,j).lcp;
                       lsbi(i,j) = nodes(i,j).lsbi;
                       gcp(i,j)  = nodes(i,j).gcp;
                       gsbi(i,j) = nodes(i,j).gsbi;  
                       p_c_l(i,j) = nodes(i,j).p_c_l;
                       p_es_l(i,j) = nodes(i,j).p_es_l;                       
                       p_c_g(i,j) = nodes(i,j).p_c_g;
                       p_es_g(i,j) = nodes(i,j).p_es_g;   
                       alpha(i,j)  = max(log10(nodes(i,j).alpha_l),1);                    
                   end
                   
                    lsbi_helper = lsbi(lsbi~=0);
                end
                figure('name','local critical points'); h(1) = imagesc(lcp); grid on;                  
                figure('name','local superbasins');     h(2) = imagesc(lsbi); grid on;  
                figure('name','gobal critical points'); h(3) = imagesc(gcp); grid on;
                figure('name','global superbasins');    h(4) = imagesc(gsbi); grid on; 
                figure('name','alpha');    h(5) = imagesc(alpha); grid on; 
                output_path = '/home/manuel/Schreibtisch/master_thesis/Thesis/imagery/';                
                figure('name','local critical probabilities');
                set(gcf,'units','normalized','outerposition',[0.25 0.2 0.5 0.71]);
                surf(X,Y,p_c_l);
                grid on; set(gca,'fontsize',20); ax = gca; ax.TickLabelInterpreter = 'latex'; 
                xlabel('position $x$ (nm)','interpreter','latex','fontsize',20);
                ylabel('position $y$ (nm)','interpreter','latex','fontsize',20);
                cb = colorbar; cb.TickLabelInterpreter = 'latex'; cb.FontSize = 20; 
                colormap(flipud(hot)); caxis([0.75 1]); xlim([1 100]); ylim([1 100]);  
                matlab2tikz([output_path,'/local_scaling_map.tex'],'width','\fwidth','height','\fheight','floatFormat','%.4f');
                figure('name','global critical probabilities');
                set(gcf,'units','normalized','outerposition',[0.25 0.2 0.5 0.71]);
                surf(X,Y,p_c_g);
                grid on; set(gca,'fontsize',20); ax = gca; ax.TickLabelInterpreter = 'latex'; 
                xlabel('position $x$ (nm)','interpreter','latex','fontsize',20);
                ylabel('position $y$ (nm)','interpreter','latex','fontsize',20);
                cb = colorbar; cb.TickLabelInterpreter = 'latex'; cb.FontSize = 20; 
                colormap jet; caxis([0.75 1]); xlim([1 100]); ylim([1 100]);                  
                matlab2tikz([output_path,'/global_scaling_map.tex'],'width','\fwidth','height','\fheight','floatFormat','%.4f');
            end            
        end   
        
        % setter 
        function nodes = set_PES(nodes,PES)
            if nargin ~= 0
                m = size(nodes,1);
                n = size(nodes,2);
                for i = 1:m
                    for j = 1:n
                        nodes(i,j).pes = PES(i,j);
                    end
                end                
            end
        end
        
   end
   
end