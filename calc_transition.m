%% Calculation of charge carrier position after hopping transition
%  input parameter: pos -> current position of charge carrier
%                   trans -> transition vector
%  output parameter: pos_t -> position of charge carrier after transition
function pos_t = calc_transition(pos,trans)

c = pp_constants;
         
if(((pos(2) + trans(2)) == 0) || ((pos(2) + trans(2)) == (nx(c) + 1)))
    % prevent hopping of holes into cathode/anode contact
    pos_t = pos;
else
    % calculate position of possible transition
    if((pos(1) + trans(1)) == 0)
        % periodic boundary condition: y down
        pos_t = [ny(c), pos(2) + trans(2)];
    elseif(pos(1) + trans(1) == (ny(c) + 1))
        % periodic boundary condition: y up
        pos_t = [1, pos(2) + trans(2)];
    else
        % inner grid-point transition
        pos_t = pos + trans;
    end  
end

end