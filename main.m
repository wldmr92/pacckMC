%% Acceleration Scheme for Charge Transport in kinetic Monte Carlo Simulations

%  Development of an accelerated kMC algorithm for the efficient treatment
%  of energetic superbasins, boundaries (amorph-crystalline) and
%  interfaces for charge transport simulations

%  This file simulates two structures (SIM1 and SIM2) which have been used
%  in Kaiser, Gößwein, Gagliardi *** SUBMITTED to Journal of Chemical Physics *** 

%% Simulation set-up
c   = pp_constants;
lut = lookup_table({'y_u';'x_r';'y_d';'x_l'},[1 0; 0 1; -1 0; 0 -1]);
% node identifier
id = reshape(1:c.x_size*c.y_size,c.y_size,c.x_size);
% lattice coordinates
x = c.r_L:c.r_L:c.x_size;
y = c.r_L:c.r_L:c.y_size;
[X,Y] = meshgrid(x,y); 

%% material and structure
mat = cell(ny(c),nx(c)); stct = mat;
structure = input('SIM1 or SIM2? [SIM1/SIM2]: ','s');
if (strcmp(structure,'SIM1') == 1)
    [mat{:,:}] = deal('P3HT'); [stct{:,:}] = deal('amorph');
elseif (strcmp(structure,'SIM2') == 1)
    [mat{:,1:idx_x(c)}] = deal('P3HT'); [stct{:,1:idx_x(c)}] = deal('amorph');    
    [mat{:,idx_x(c)+1:2*idx_x(c)}] = deal('P3HT');   [stct{:,idx_x(c)+1:2*idx_x(c)}] = deal('crystalline');   
    [mat{:,2*idx_x(c)+1:3*idx_x(c)}] = deal('P3HT'); [stct{:,2*idx_x(c)+1:3*idx_x(c)}] = deal('amorph');    
    [mat{:,3*idx_x(c)+1:4*idx_x(c)}] = deal('P3HT'); [stct{:,3*idx_x(c)+1:4*idx_x(c)}] = deal('crystalline');  
    [mat{:,4*idx_x(c)+1:end}] = deal('P3HT');        [stct{:,4*idx_x(c)+1:end}] = deal('amorph'); 
else
    warning('Unknown structure!!!'); return;
end
area_a = (strcmp(stct,'amorph') == 1);
area_c = (strcmp(stct,'crystalline') == 1);

%% energetic disorder (diagonal and off-diagonal)
% diagonal
sigma_a = [0; 50e-3; 100e-3];
sigma_c = sigma_a./2;
sigma = [sigma_a, sigma_c];
[n_sigma,~] = size(sigma);
% off-diagonal
sigma_od = ((2.*c.gamma.*c.r_L)/sqrt(2)).*[0; 0.5; 1];
n_sigma_od = numel(sigma_od);
% prefactor
a = c.a_0.*area_a + c.a_1*area_c;
%% bias voltage
V_bias = transpose(logspace(-1,2,10));
n_bias = numel(V_bias);
% eletrical field strength (V/cm)
E = (V_bias./c.x_size).*1e7;
%% output parameters: time of flight variables:
t_tof  = zeros(c.n_sim,n_bias);
v_av   = zeros(n_bias,n_sigma_od);
mu_av  = zeros(n_bias,n_sigma_od);
%% output parameters: performance measurements
t_sim  = zeros(n_sigma_od,n_sigma);
t_set_up = zeros(n_sigma_od,n_sigma);
t_sim_kMC  = zeros(c.n_sim,n_bias);
n_sim_kMC  = zeros(c.n_sim,n_bias);
%% Generate simulation set-up
% diagonal energetic disorder
prompt = input('Generate new diagonal disorder? [Y/N]: ','s');
for l = 1:n_sigma
    if (strcmp(prompt,'Y') == 1)
        if exist(['./sim_set_up/',structure,'/E_i_HOMO_',num2str(sigma(l)),'.txt'], 'file')
            delete(['./sim_set_up/',structure,'/E_i_HOMO_',num2str(sigma(l)),'.txt']);             
        end     
        rng('shuffle');
        E_i_HOMO = c.E_0_HOMO.*area_a + sigma(l,1).*randn(ny(c),nx(c)).*area_a + c.E_0_HOMO.*area_c + sigma(l,2).*randn(ny(c),nx(c)).*area_c;
        save(['./sim_set_up/',structure,'/E_i_HOMO_',num2str(sigma(l)),'.txt'],'E_i_HOMO','-ascii');
    elseif(strcmp(prompt,'N') == 1)
        disp('No new diagonal disorder generated!'); break;
    else
        warning('Unknown input!!!'); return;
    end
end
% off-diagonal energetic disorder
 prompt = input('Generate new off-diagonal disorder? [Y/N]: ','s');
 for l = 1:n_sigma_od
     if (strcmp(prompt,'Y') == 1)
         if exist(['./sim_set_up/',structure,'/Gamma_sigma_od_',num2str(sigma_od(l)),'.txt'], 'file')
            delete(['./sim_set_up/',structure,'/Gamma_sigma_od_',num2str(sigma_od(l)),'.txt']);           
         end 
         rng('shuffle');           
         Gamma = (c.gamma.*c.r_L + sigma_od(l).*randn(ny(c),nx(c))).*area_a  + c.gamma.*c.r_L.*area_c;
         save(['./sim_set_up/',structure,'/Gamma_sigma_od_',num2str(sigma_od(l)),'.txt'],'Gamma','-ascii');
     elseif(strcmp(prompt,'N') == 1)
         disp('No new off-diagonal disorder generated!'); break;
     else
        warning('Unknown input!!!'); return;
     end
 end
% simulation type: accelerated/non-accelerated
prompt = input('Apply Scaling Algorithm? [Y/N]: ','s');
if (strcmp(prompt,'Y') == 1)
    scaling = 'scaling_on';
elseif(strcmp(prompt,'N') == 1)
    scaling = 'scaling_off';
else
    warning('Unknown input!!!'); return;
end
% generate file  for output data
if (strcmp(structure,'SIM1') == 1 )
    if (strcmp(scaling,'scaling_on') == 1 )
        FE = '_as_a.txt';
    elseif (strcmp(scaling,'scaling_off') == 1 )
        FE = '_nas_a.txt';        
    end
elseif (strcmp(structure,'SIM2') == 1)
    if (strcmp(scaling,'scaling_on') == 1 )
        FE = '_as_c.txt';        
    elseif (strcmp(scaling,'scaling_off') == 1 )
        FE = '_nas_c.txt';        
    end    
end
%% input parameters: Acceleration Scheme 
input('Specify input parameters for acceleration scheme: [Press ENTER!] ');
% minimal critical process propability for defining superbasins
xi = input('Minimal critical process propability for defining superbasins: xi = ');
% error measure
delta = input('Error Measure: delta = ');  
% growth rate of scaling factor in case of multiple scaling
kappa = input('Growth Rate of Scaling Factor (alpha^(m/kappa)): kappa = ');  
% variables for storing superbasin times
t_sbi_l = cell(c.n_sim,1); t_sbi_g = cell(c.n_sim,1); 
% variables for storing superbasin counter values
n_sbi_l = cell(c.n_sim,1); n_sbi_g = cell(c.n_sim,1); 
% variables for storing superbasin multiple scaling factor
m_sbi_l = cell(c.n_sim,1); m_sbi_g = cell(c.n_sim,1); 
%% Generate datapath for saving results
datapath = ['./results/',structure,'/xi_',num2str(xi),'_delta_',num2str(delta),'_kappa_',num2str(kappa)]; 

if ~exist(datapath, 'dir')
    mkdir(datapath);
    % generate folder for output data
    output = {'t_tof'; 'v_av'; 'mu_av'; 't_sim'; 't_sbi'};
    for l = 1:numel(output)
        mkdir([datapath,'/',output{l}]);
    end
end
 %% Simulation Sweeps
 % on-diagonal disorder sweep
 for j = n_sigma:-1:1
     % off-diagonal disorder sweep
     for l = n_sigma_od:-1:1
         % first timer: duration of overall simulation set-up including
         % scaling procedures
         tic
         % load diagonal energetic disorder
         E_i_HOMO = load(['./sim_set_up/',structure,'/E_i_HOMO_',num2str(sigma(j,1)),'.txt']);
         % load off-diagonal energetic disorder
         Gamma = load(['./sim_set_up/',structure,'/Gamma_sigma_od_',num2str(sigma_od(l)),'.txt']);
         % potential energy surface (at V_bias = 0)
         PES = E_i_HOMO;
         % construct lattice of nodes
         lattice = nodes(id,mat,stct,X,Y,PES,a,Gamma);
         % global superbasin reference energy barrier
         if (strcmp(structure,'SIM1') == 1)
            a_g = c.a_0.*exp(-2.*(c.gamma - 2.*sigma_od(l)).*c.r_L);
         elseif (strcmp(structure,'SIM2') == 1)
            a_g = sqrt(c.a_0*c.a_1).*exp(-2.*c.gamma.*c.r_L);             
         end       
         E_a_ref_g = calc_equivalent_activation_energy(a_g);
         % calculate static rate constants and global critical points
         [lattice,E_a,GAMMA] = calc_static_rates(lattice,lut,E_a_ref_g,'scaling_on');
         % compute local critical points
         lattice = calc_critical_points(lattice,lut,xi,'hole');
         % group suberbasins
         [lattice,idx_g,idx_l] = group_superbasins(lut,lattice); 
         % average occupation time per superbasin
         t_sbi_l_av  = zeros(idx_l,n_bias); t_sbi_g_av = zeros(idx_g,n_bias);
         t_sbi_l_std = zeros(idx_l,n_bias); t_sbi_g_std = zeros(idx_g,n_bias);
         % average counter value per superbasin
         n_sbi_l_av = zeros(idx_l,n_bias); n_sbi_g_av = zeros(idx_g,n_bias);
         n_sbi_l_std = zeros(idx_l,n_bias); n_sbi_g_std = zeros(idx_g,n_bias);
         % average superbasin multiple scaling factor
         m_sbi_l_av = zeros(idx_l,n_bias); m_sbi_g_av = zeros(idx_g,n_bias);
         m_sbi_l_std = zeros(idx_l,n_bias); m_sbi_g_std = zeros(idx_g,n_bias);         
         % calc scaling factors
         lattice = calc_scaling_factors(lattice,idx_g,idx_l,xi,delta,kappa,'hole'); 
         display(lattice,X,Y);
         % number of critical sightings
         
         N_fl = zeros(idx_l,2); N_fg = zeros(idx_g,2);
         % generate charge carrier object
         cc = particles({'hole'},1,N_fl,N_fg);
         t_set_up(l,j) = toc; save([datapath,'/t_sim/t_set_up',FE],'t_set_up','-ascii');
         % second timer: overall simulation time -> bias sweep for fixed on/off diagonal disorder
         tic
         for k = 1:n_bias
             % potential energy surface
             E_i_bias = (X./c.x_size).*V_bias(k);
             PES = E_i_bias + E_i_HOMO;
             % set potential energy surface
             lattice = set_PES(lattice,PES);
             % calculate static rate constants and global critical points
             [lattice,~] = calc_static_rates(lattice,lut,E_a_ref_g,'scaling_off');
             parfor i = 1:c.n_sim
                 disp('Start kMC run'); disp(i)
                 %% Calculation of rates, Monte Carlo Step, System Update; Scaling, Termination Condition
                 [t_tof(i,k),t_sim_kMC(i,k),n_sim_kMC(i,k),t_sbi_l{i},t_sbi_g{i},n_sbi_l{i},n_sbi_g{i},m_sbi_l{i},m_sbi_g{i}] = apply_kMC(lattice,lut,cc,idx_g,idx_l,kappa,scaling,delta);
             end
             % average occupation time per superbasin
             t_sbi_l_av(:,k) = mean(transpose([t_sbi_l{:}])); t_sbi_l_std(:,k) = std(transpose([t_sbi_l{:}]));
             t_sbi_g_av(:,k) = mean(transpose([t_sbi_g{:}])); t_sbi_g_std(:,k) = std(transpose([t_sbi_g{:}]));
             % average counter value per superbasin
             n_sbi_l_av(:,k) = mean(transpose([n_sbi_l{:}])); n_sbi_l_std(:,k) = std(transpose([n_sbi_l{:}]));
             n_sbi_g_av(:,k) = mean(transpose([n_sbi_g{:}])); n_sbi_g_std(:,k) = std(transpose([n_sbi_g{:}])); 
             % variables for storing superbasin multiple scaling factor
             m_sbi_l_av(:,k) = mean(transpose([m_sbi_l{:}])); m_sbi_l_std(:,k) = std(transpose([m_sbi_l{:}]));
             m_sbi_g_av(:,k) = mean(transpose([m_sbi_g{:}])); m_sbi_g_std(:,k) = std(transpose([m_sbi_g{:}]));                         
             %% Post-Processing
             % average drift velocity
             v_av(k,l) = exp(transpose(mean(log(1./t_tof(:,k))))).*c.x_size.*1e-7;
             % average mobility
             mu_av(k,l) = v_av(k,l)./E(k);
             % save results
             save([datapath,'/t_tof/t_tof_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'t_tof','-ascii');
             save([datapath,'/t_sim/t_sim_kMC_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'t_sim_kMC','-ascii');
             save([datapath,'/t_sim/n_sim_kMC_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'n_sim_kMC','-ascii');
             save([datapath,'/t_sbi/t_sbi_l_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'t_sbi_l_av','-ascii');
             save([datapath,'/t_sbi/t_sbi_g_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'t_sbi_g_av','-ascii');
             save([datapath,'/t_sbi/n_sbi_l_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'n_sbi_l_av','-ascii');
             save([datapath,'/t_sbi/n_sbi_g_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'n_sbi_g_av','-ascii');  
             save([datapath,'/t_sbi/m_sbi_l_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'m_sbi_l_av','-ascii');
             save([datapath,'/t_sbi/m_sbi_g_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'m_sbi_g_av','-ascii');               
             save([datapath,'/t_sbi/t_sbi_l_std_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'t_sbi_l_std','-ascii');
             save([datapath,'/t_sbi/t_sbi_g_std_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'t_sbi_g_std','-ascii');
             save([datapath,'/t_sbi/n_sbi_l_std_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'n_sbi_l_std','-ascii');
             save([datapath,'/t_sbi/n_sbi_g_std_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'n_sbi_g_std','-ascii');   
             save([datapath,'/t_sbi/m_sbi_l_std_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'m_sbi_l_std','-ascii');
             save([datapath,'/t_sbi/m_sbi_g_std_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'m_sbi_g_std','-ascii');                
             save([datapath,'/v_av/v_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'v_av','-ascii');
             save([datapath,'/mu_av/mu_av_sigma_',num2str(sigma(j)),'_SIGMA_',num2str(sigma_od(l)),FE],'mu_av','-ascii');
         end
         t_sim(j,l) = toc; save([datapath,'/t_sim/t_sim',FE],'t_sim','-ascii');
     end
     % reset output variables: time of flight
     v_av   = zeros(n_bias,n_sigma_od);
     mu_av  = zeros(n_bias,n_sigma_od);
 end