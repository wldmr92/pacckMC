%% Kinetic Monte Carlo Method
%  input parameter: lattice -> lattice object
%                   lut -> lookup table object for rate catalog
%                   cc -> charge carrier object
%                   idx_g -> number of global superbasins
%                   idx_l -> number of local superbasins
%                   kappa -> growth rate of scaling factor
%                   scaling -> string for on/off switching of scaling algorithm
%                   delta -> uncertainty
%  output parameter: t_tof -> time of flight of kMC simulation
%                    t_sim_kMC -> simulation time of one kMC look
%                    n_sim_kMC -> number of kMC loop iterations
%                    t_lsbi -> occupation time of local superbasins
%                    t_gsbi -> occupation time of global superbasins 
%                    n_lsbi -> counter for local superbasins
%                    n_gsbi -> counter for global superbasins 
function [t_tof,t_sim_kMC,n_sim_kMC,t_lsbi,t_gsbi,n_lsbi,n_gsbi,m_lsbi,m_gsbi] = apply_kMC(lattice,lut,cc,idx_g,idx_l,kappa,scaling,delta)

c = pp_constants;  
% number of possible transitions
n_trans = numel(lut);
% occupation time superbasins
t_lsbi = zeros(idx_l,1);
t_gsbi = zeros(idx_g,1);
% counter superbasins
n_lsbi = zeros(idx_l,1);
n_gsbi = zeros(idx_g,1);
% counter superbasin multiple scaling factor
m_lsbi = zeros(idx_l,1);
m_gsbi = zeros(idx_g,1);
% counter for kMC-loop interations
n_sim_kMC = 0;

% second timer: simulation time per kMC loop
tic
 %% Initialization of kMC-Loop
 rng('shuffle');
 % occupation of lattice
 cc.pos = [randi([1 ny(c)]), 1]; pos = cc.pos;
 % reset current transition and time
 mu = 0; t = 0;     
 % store transitions
 trans = {lut.trans};
 if(strcmp(scaling,'scaling_on') == 1)
     %% kMC-Loop
     while(1)
         % update counter for kMC-loop interations
         n_sim_kMC = n_sim_kMC + 1;
         % get static rates **** CALCULATION OF RATES **** 
         src = lattice(pos(1),pos(2)).src_f(:);
         %% apply scaling
         if (cc.alpha_l == 1 && cc.alpha_g == 1)
             for k = 1:n_trans
                 lut(k).rate = src(k);
                 % calculate partial sums
                 lut(k).psum = sum([lut(1:k).rate]);
             end
         % prioritize local scaling factor
         elseif cc.alpha_l > 1
             % get local critical branches
             lcb = lattice(pos(1),pos(2)).lcb(:);
             % fill rate cataloge
             for k = 1:n_trans
                 if( lcb(k) == 2 )
                     lut(k).rate = src(k)/cc.alpha_l;
                 elseif( lcb(k) == 1 )    
                     % get position of possible transition
                     pos_t = calc_transition(pos,trans{k});
                     if((cc.lsbi(1) == 0 && lattice(pos_t(1),pos_t(2)).lsbi == cc.lsb) || (lattice(pos_t(1),pos_t(2)).lsbi == 0 && cc.lsbi(1) > 0))
                        lut(k).rate = src(k);
                     else
                        lut(k).rate = src(k);  
                     end
                 else
                     lut(k).rate = src(k);
                 end
                 % calculate partial sums
                 lut(k).psum = sum([lut(1:k).rate]);
             end
         elseif cc.alpha_g > 1
             % get global critical branches
             gcb = lattice(pos(1),pos(2)).gcb(:);
             % fill rate cataloge
             for k = 1:n_trans
                 if( gcb(k) == 2 )
                     lut(k).rate = src(k)/cc.alpha_g;
                 elseif( gcb(k) == 1 )
                     % get position of possible transition
                     pos_t = calc_transition(pos,trans{k});
                     if((cc.gsbi(1) == 0 && lattice(pos_t(1),pos_t(2)).gsbi == cc.gsb) || (lattice(pos_t(1),pos_t(2)).gsbi == 0 && cc.gsbi(1) > 0))
                        lut(k).rate = src(k);
                     else
                        lut(k).rate = src(k);  
                     end
                 else
                     lut(k).rate = src(k);
                 end
                 % calculate partial sums
                 lut(k).psum = sum([lut(1:k).rate]);
             end
         end
         
         %% Monte Carlo Step
         % draw two uniform random numbers
         r_1 = rand();
         r_2 = rand();
         % get total rate
         a_tot = lut(end).psum;
         % generate array of partial sums
         a_psum = [lut(:).psum];
         % multiply the total rate by r_1
         a_mu = r_1.*a_tot;
         % determine transition
         mu = sum(a_psum < a_mu) + 1;
         % calculate time step
         tau = -log(r_2)./a_tot;
         
         %% System Update
         % update system time
         t = t + tau;
         % execute transition
         cc.pos = calc_transition(cc.pos,lut(mu).trans);
         pos = cc.pos;
         
         %% global superbasins
         if idx_g > 0
             % store the last two global superbasin indices
             cc.gsbi(2:3) = cc.gsbi(1:2);
             % get current global superbasin index
             cc.gsbi(1) = lattice(pos(1),pos(2)).gsbi;
             gsbi = cc.gsbi(1);
             % update global superbasin counter
             if gsbi ~= 0
                 n_gsbi(gsbi) = n_gsbi(gsbi) + 1;
             end
             % update global superbasin time
             if cc.gsbi(2) ~= 0
                 t_gsbi(cc.gsbi(2)) = t_gsbi(cc.gsbi(2)) + tau;
             end
             % check whether a new global superbasin was entered
             if(cc.gsbi(1) > 0 && ( (cc.gsbi(3) > 0 && cc.gsbi(1) ~= cc.gsbi(3)) || (cc.gsbi(2) == 0 && cc.gsbi(3) == 0)))
                 % set current global superbasin
                 cc.gsb = gsbi;
             end
             % increment critical counter and set scaling factor
             if(cc.gsbi(1) > 0)
                 % increment counter for sightings in local superbasin
                 cc.N_fg(gsbi,1) = cc.N_fg(gsbi,1) + 1;
                 % get number of critical sightings in local superbasin
                 N_fg = update_cc_N_f(delta,lattice(pos(1),pos(2)).alpha_g, cc.N_fg(gsbi,2), kappa);
                 % check if counter has reached number of critical sightings
                 if (cc.N_fg(gsbi,1) > N_fg)%(mod(cc.N_fg(gsbi,1),N_fg) == 0)
                     % increment how often the number of critical sightings has
                     % already been reached in the current superbasin
                     cc.N_fg(gsbi,2) = cc.N_fg(gsbi,2) + 1;
					 cc.N_fg(gsbi,1) = 0;
                     m_gsbi(gsbi) = m_gsbi(gsbi) + 1;
                     % get scaling factor
                     cc.alpha_g = update_cc_alpha(lattice(pos(1),pos(2)).alpha_g, cc.N_fg(gsbi,2), kappa);
                 end
             end
             % reset global scaling factor
             if (cc.gsbi(1) == 0 && cc.gsbi(2) == 0 && cc.alpha_g > 1)
                 for k = 1:n_trans
                     % check wether all nearest neighbors are non-critical nodes
                     pos_t = calc_transition(cc.pos,lut(k).trans);
                     if(lattice(pos_t(1),pos_t(2)).gsbi == 0)
                         if k == n_trans
                             cc.alpha_g = 1;
                             if (cc.gsbi(3) > 0)
                                 cc.N_fg(cc.gsbi(3),1) = 0;
                                 cc.N_fg(cc.gsbi(3),2) = 0;
                             end
                         end
                     else
                         break;
                     end
                 end
             end
         end
         
         %% local superbasins
         if idx_l > 0
             % store the last two local superbasin indices
             cc.lsbi(2:3) = cc.lsbi(1:2);
             % get current local superbasin index
             cc.lsbi(1) = lattice(pos(1),pos(2)).lsbi;
             lsbi = cc.lsbi(1);
             % update local superbasin counter
             if lsbi ~= 0
                 n_lsbi(lsbi) = n_lsbi(lsbi) + 1;
             end
             % update local superbasin time
             if cc.lsbi(2) ~= 0
                 t_lsbi(cc.lsbi(2)) = t_lsbi(cc.lsbi(2)) + tau;
             end
             % check wether a new local superbasin was entered
             if(cc.lsbi(1) > 0 && ( (cc.lsbi(3) > 0 && cc.lsbi(1) ~= cc.lsbi(3)) || (cc.lsbi(2) == 0 && cc.lsbi(3) == 0)))
                 % set current local superbasin
                 cc.lsb = lsbi;         
             end
             % increment critical counter and set scaling factor
             if(cc.lsbi(1) > 0)
                 % increment counter for sightings in local superbasin
                 cc.N_fl(lsbi,1) = cc.N_fl(lsbi,1) + 1;
                 % get number of critical sightings in local superbasin
                 N_fl = update_cc_N_f(delta, lattice(pos(1),pos(2)).alpha_l, cc.N_fl(lsbi,2), kappa);
                 % check if counter has reached number of critical sightings
                 if (cc.N_fl(lsbi,1) > N_fl)
                     % increment how often the number of critical sightings has
                     % already been reached in the current superbasin
					 cc.N_fl(lsbi,1) = 0; 
                     cc.N_fl(lsbi,2) = cc.N_fl(lsbi,2) + 1; 
                     m_lsbi(lsbi) = m_lsbi(lsbi) + 1;
                     % get scaling factor
                     cc.alpha_l = update_cc_alpha(lattice(pos(1),pos(2)).alpha_l, cc.N_fl(lsbi,2), kappa);
                 end
             end
             % reset local scaling factor
             if (cc.lsbi(1) == 0 && cc.lsbi(2) == 0 && cc.alpha_l > 1)
                 for k = 1:n_trans
                     % check wether all nearest neighbors are non-critical nodes
                     pos_t = calc_transition(cc.pos,lut(k).trans);
                     if(lattice(pos_t(1),pos_t(2)).lsbi == 0)
                         if k == n_trans
                             cc.alpha_l = 1;
                             if (cc.lsbi(3) > 0)
                                 cc.N_fl(cc.lsbi(3),1) = 0;
                                 cc.N_fl(cc.lsbi(3),2) = 0;
                             end
                         end
                     else
                         break;
                     end
                 end
             end
         end
         
         %% Check termination condition: hole reached cathode
         if(cc.pos(2) == nx(c))
             % get time of flight 
             t_tof = t; fprintf('t_tof = %4.8f s\n',t_tof);
             % second timer: save simulation time per kMC loop
             t_sim_kMC = toc; fprintf('t_sim_kMC = %4.4f s\n',t_sim_kMC);        
             break;
         end           
         
     end
 elseif(strcmp(scaling,'scaling_off') == 1)
     %% kMC-Loop
     while(1)
         % update counter for kMC-loop interations
         n_sim_kMC = n_sim_kMC + 1;
         % get static rates
         src = lattice(pos(1),pos(2)).src_f(:);
         % fill rate cataloge
         for k = 1:n_trans
             lut(k).rate = src(k);
             % calculate partial sums
             lut(k).psum = sum([lut(1:k).rate]);
         end
         
         %% Monte Carlo Step
         % draw two uniform random numbers
         r_1 = rand();
         r_2 = rand();
         % get total rate
         a_tot = lut(end).psum;
         % generate array of partial sums
         a_psum = [lut(:).psum];
         % multiply the total rate by r_1
         a_mu = r_1.*a_tot;
         % determine transition
         mu = sum(a_psum < a_mu) + 1;
         % calculate time step
         tau = -log(r_2)./a_tot;
         
         %% System Update
         % update system time
         t = t + tau;
         % execute transition
         cc.pos = calc_transition(cc.pos,lut(mu).trans);
         pos = cc.pos;

         %% global superbasins
         if idx_g > 0
             % store the last two global superbasin indices
             cc.gsbi(2:3) = cc.gsbi(1:2);
             % get current global superbasin index
             cc.gsbi(1) = lattice(pos(1),pos(2)).gsbi;
             gsbi = cc.gsbi(1);
             % update global superbasin counter
             if gsbi ~= 0
                 n_gsbi(gsbi) = n_gsbi(gsbi) + 1;
             end
             % update global superbasin time
             if cc.gsbi(2) ~= 0
                 t_gsbi(cc.gsbi(2)) = t_gsbi(cc.gsbi(2)) + tau;
             end
         end
         
         %% local superbasins
         if idx_l > 0
             % store the last two local superbasin indices
             cc.lsbi(2:3) = cc.lsbi(1:2);
             % get current local superbasin index
             cc.lsbi(1) = lattice(pos(1),pos(2)).lsbi;
             lsbi = cc.lsbi(1);
             % update local superbasin counter
             if lsbi ~= 0
                 n_lsbi(lsbi) = n_lsbi(lsbi) + 1;
             end
             % update local superbasin time
             if cc.lsbi(2) ~= 0
                 t_lsbi(cc.lsbi(2)) = t_lsbi(cc.lsbi(2)) + tau;
             end
         end
         
         %% Check termination condition: hole reached cathode
         if(cc.pos(2) == nx(c))
             % get time of flight 
             t_tof = t; fprintf('t_tof = %4.8f s\n',t_tof);
             % second timer: save simulation time per kMC loop
             t_sim_kMC = toc; fprintf('t_sim_kMC = %4.4f s\n',t_sim_kMC);        
             break;
         end           
     end
 else
     warning('Unknown input parameter!!!'); return;
 end
     
end