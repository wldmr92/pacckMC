%% Convert rate constant in equivalent activation energy of Arrhenius equation
%  input parameter: k_ij -> rate constant 
%  output parameter: E_a -> equivalent activation energy
function E_a = calc_equivalent_activation_energy(k_ij)

c = pp_constants;
% thermal energy
E_th = c.k_B.*c.T;
% prefactor reference Arrhenius equation
nu = E_th/c.h;
% equivalent energy barrier of Arrhenius equation
E_a = E_th*log(k_ij./nu);

end