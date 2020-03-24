%% Calculate partition function of Boltzmann distribution for quasi-equilibrated superbasin
%  input parameters: F_i -> free energies 
%  ouput parameters: Z -> partition function
function Z = calc_partition_function(F_i)

c = pp_constants;
% thermal energy
E_th = c.k_B.*c.T;
% partition function
Z = sum(exp(-F_i./E_th));

end