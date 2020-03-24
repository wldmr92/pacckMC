function  R_ij = calc_distance(pos,pos_t)

c = pp_constants;

Delta_y1 = abs(pos_t(1) - pos(1));
Delta_y2 = c.y_size - Delta_y1;

Delta_x = abs(pos(2) - pos_t(2));

if(Delta_y1 <= Delta_y2)
    R_ij = sqrt( Delta_x.^2 + Delta_y1.^2 );
else
    R_ij = sqrt( Delta_x.^2 + Delta_y2.^2 );    
end

end