function img_new = downsample_image(img, N_new)

    % img should be a N by N square
    N = size(img, 1);
    img_new = zeros(N_new, N_new);
    
    unit_new = N/N_new;
    
    for (i=1:N)
        for (j=1:N)
            i_left = ceil((i-1)/unit_new);
            i_right = ceil(i/unit_new);
            j_up = ceil((j-1)/unit_new);
            j_down = ceil(j/unit_new);
            
            if ((i_left==i_right) && (j_up==j_down))
                if ((i_left~=0) && (j_up~=0))
                img_new(i_left,j_up) = img_new(i_left,j_up) + img(i,j);
                end
            elseif ((i_left==i_right) && (j_up~=j_down))
                j_up_part = j_up*unit_new-(j-1);
                if ((i_left~=0) && (j_up~=0))
                img_new(i_left,j_up) = img_new(i_left,j_up) + j_up_part * img(i,j);
                end
                if ((i_left~=0) && (j_down~=0))
                img_new(i_left,j_down) = img_new(i_left,j_down) + (1-j_up_part) * img(i,j);
                end
            else if ((i_left~=i_right) && (j_up==j_down))
                i_left_part = i_left*unit_new-(i-1);
                if ((i_left~=0) && (j_up~=0))
                img_new(i_left,j_up) = img_new(i_left,j_up) + i_left_part * img(i,j);
                end
                if ((i_right~=0) && (j_up~=0))
                img_new(i_right,j_up) = img_new(i_right,j_up) + (1-i_left_part) * img(i,j);
                end
            else
                i_left_part = i_left*unit_new-(i-1);
                j_up_part = j_up*unit_new-(j-1);
                
                if ((i_left~=0) && (j_up~=0))
                img_new(i_left, j_up) = img_new(i_left,j_up) + i_left_part*j_up_part * img(i,j);
                end
                if ((i_left~=0) && (j_down~=0))
                img_new(i_left, j_down) = img_new(i_left,j_down) + i_left_part*(1-j_up_part) * img(i,j);
                end
                if ((i_right~=0) && (j_up~=0))
                img_new(i_right, j_up) = img_new(i_right,j_up) + (1-i_left_part)*j_up_part * img(i,j);
                end
                if ((i_right~=0) && (j_down~=0))
                img_new(i_right, j_down) = img_new(i_right,j_down) + (1-i_left_part)*(1-j_up_part) * img(i,j);
                end
            end
        end
    end

end
