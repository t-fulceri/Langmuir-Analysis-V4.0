obj = plasma_ST_profiles;
x_pos = 11;
y_pos = 9;
nr = 1;

cond_x = (x_pos >= obj.x_begin && x_pos <= obj.x_end);
cond_y = (y_pos >= obj.y_begin && y_pos <= obj.y_end);
cond_rep = (nr >= 1 && nr <= obj.Nr);

if(cond_x && cond_y && cond_rep)
    [~,nx] = min(abs(obj.x_axis - x_pos));
    [~,ny] = min(abs(obj.y_axis - y_pos));
    pp = obj.pp_ca{nx,ny,nr};
    return
else
    fprintf('Error: requested position or repetition index is out of range.\n')
end