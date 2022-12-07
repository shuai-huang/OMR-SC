function X = At_op_3d(y, s_vect, mat_sz)   % assume square region with a odd-number size
    X_tmp = zeros(mat_sz);
    for (i=1:length(y))
        s_vect_tmp = s_vect{i};
        X_tmp(s_vect_tmp) = X_tmp(s_vect_tmp) + y(i);
    end
    X = X_tmp;   % comples images do not need fftshift??? 
end

