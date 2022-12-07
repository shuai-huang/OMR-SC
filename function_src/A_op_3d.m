% construct linear measurement operator and its transpose operator
function y = A_op_3d(X, s_vect)
    y = zeros(length(s_vect),1);
    for (i=1:length(s_vect))
        s_vect_tmp = s_vect{i};
        y(i) = sum(X(s_vect_tmp));
    end
end

