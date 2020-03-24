%% Calculate index of reverse transitions
%  input parameters: k -> transition index of forward transition (i -> j)
%  output parameters: kp -> transitions index of backward transtion (j -> i)
function kp = calc_reverse_idx(k)

if k < 3
    kp = k+2;
else
    kp = k-2;
end

end