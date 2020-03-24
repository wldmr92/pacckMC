classdef pp_constants
%% physical constants and programmatic parameters
    properties(Constant = true)
            e = 1.6021766208e-19;       % elementary charge (As)
            c = 299792458e9;            % vacuum speed of light (nm/s)            
            k_B = 8.6173303e-5;         % Boltzmann constant (eV/K)
            h = 4.135667662e-15;        % Planks constant (eV.s)
            E_0_HOMO = -5;              % mean molecular orbital energy HOMO,P3HT (eV)
            E_0_LUMO = -3.8;            % mean molecular orbital energy LUMO,PCBM (eV)
            eps_r = 5;                  % relative permittivity (1)
            T = 300;                    % absolute temperature (K)
            a_0 = 1e10;                 % attempt-to-hop frequency amorph (1/s)
            a_1 = 1e13;                 % attempt-to-hop frequency crystalline (1/s)
            gamma = 2;                  % localization constant (1/nm)
            r_L = 1;                    % average lattice spacing (nm)
            x_size = 100;               % grid size in x-direction (nm)
            y_size = 100;               % grid size in y-direction (nm)
            n_sim = 1e2;                % number of simulations (1)
    end
    methods
        function r = mu_0(c)            % vacuum permeability (Vs/A.nm)
            r = 4.*pi.*1e-16;
        end       
        function r = eps_0(c)           % vacuum permittivity (As/V.nm)
            r = 1./([mu_0(c)].*[(c.c).^2]);
        end    
        function r = nx(c)              % number of lattice points in x-direction (1)
            r = ceil([c.x_size]./[c.r_L]);
        end
        function r = ny(c)              % number of lattice points in y-direction (1)
            r = ceil([c.y_size]./[c.r_L]);
        end    
        function r = idx_x(c)           % one-fifth of lattice points in x-direction (1)
            r = ceil(nx(c)/5);
        end
        function r = idx_y(c)           % one-fifth of lattice points in y-direction (1)
            r = ceil(ny(c)/5);
        end        
    end
end