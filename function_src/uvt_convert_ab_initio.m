function foo = uvt_noisy_convert_ab_initio(density_all, N_next, N_half_next, downsample_unit, low_res_file, ab_init_file)

    % density_all:      the sum of all the density values, can be estimated from the spatial radial features
    % N_next:           an odd number, the size of the upsampling Cartesian graid is N_next by N_next by N_next
    % N_half_next:      half of N_next, i.e. (N_half-1)/2
    % downsample_unit:  the size of the downsampling block, starting from 3 by 3, 5 by 5, 7 by 7, etc. We downsample by summing up all the pixels within the downsampling block
    % low_res_file:     the recovered low-res density by OMR-SC via downsampling, to be used as the ab initio model
    % ab_init_file:     the upsampled density from low_res_file


    % N is the size of the sampling grid that encompasses the downsampled density
    N = 2*floor((N_half_next-(downsample_unit-1)/2)/downsample_unit)+1;
    N_half = (N-1)/2;

    den_3d_rec = repmat(0,[N N N]);
    den_rec = dlmread(strcat(low_res_file));

    for (i=1:size(den_rec,1))
        den_3d_tmp = den_rec(i,:);
        den_3d_rec(den_3d_tmp(1)+N_half+1,den_3d_tmp(2)+N_half+1,den_3d_tmp(3)+N_half+1) = den_3d_tmp(4);
    end

    x_ori_pre = []; 
    weight_ori_pre = []; 
    weight_thd = 1e-12;
    for (i=1:N)
        for (j=1:N)
            for (k=1:N)
                if (den_3d_rec(i,j,k)>weight_thd) 
                    x_ori_pre = [x_ori_pre; i-N_half-1 j-N_half-1 k-N_half-1];
                    weight_ori_pre = [weight_ori_pre; den_3d_rec(i,j,k)];
                end
            end
        end
    end

    %%%%%%%%%%%%
    x_ori_pre = x_ori_pre * N_next/N;
    s_sd = 0.8660254 * N_next/N;

    den_3d_next = zeros(N_next,N_next,N_next);

    s_all = x_ori_pre;

    for (s_idx=1:size(s_all,1))
        s_tmp = s_all(s_idx,:);
        weight_tmp = weight_ori_pre(s_idx);

        den_tmp_1=normpdf( (round(s_tmp(1)-3*s_sd)-1):(round(s_tmp(1)+3*s_sd)+1) , s_tmp(1), s_sd);
        den_tmp_2=normpdf( (round(s_tmp(2)-3*s_sd)-1):(round(s_tmp(2)+3*s_sd)+1) , s_tmp(2), s_sd);
        den_tmp_3=normpdf( (round(s_tmp(3)-3*s_sd)-1):(round(s_tmp(3)+3*s_sd)+1) , s_tmp(3), s_sd);

        idx_i=0;
        for (i=(round(s_tmp(1)-3*s_sd)-1):(round(s_tmp(1)+3*s_sd)+1))
            idx_i = idx_i+1;
            idx_j = 0;
            for (j=(round(s_tmp(2)-3*s_sd)-1):(round(s_tmp(2)+3*s_sd)+1))
                idx_j = idx_j+1;
                idx_k = 0;
                for (k=(round(s_tmp(3)-3*s_sd)-1):(round(s_tmp(3)+3*s_sd)+1))
                    idx_k = idx_k+1;
                    den_tmp = weight_tmp*den_tmp_1(idx_i)*den_tmp_2(idx_j)*den_tmp_3(idx_k);
                    if (i+N_half_next+1>=1)&&(i+N_half_next+1<=N_next)
                        if (j+N_half_next+1>=1)&&(j+N_half_next+1<=N_next)
                            if (k+N_half_next+1>=1)&&(k+N_half_next+1<=N_next)
                                den_3d_next(i+N_half_next+1,j+N_half_next+1,k+N_half_next+1) = den_3d_next(i+N_half_next+1,j+N_half_next+1,k+N_half_next+1) + den_tmp;
                            end
                        end
                    end
                end
            end
        end
    end



    den_3d_seq = den_3d_next;
    den_3d_out = zeros(length(den_3d_seq(den_3d_seq~=0)), 4);
    idx=1;
    for (i=1:N_next)
        for (j=1:N_next)
            for (k=1:N_next)
                if (den_3d_next(i,j,k)~=0)
                    den_3d_out(idx,:) = [i-N_half_next-1 j-N_half_next-1 k-N_half_next-1 den_3d_next(i,j,k)];
                    idx = idx+1;
                end
            end
        end
    end

    den_3d_out(:,4) = den_3d_out(:,4)/sum(den_3d_out(:,4))*density_all;
    dlmwrite(ab_init_file, den_3d_out, 'delimiter', ' ', 'precision', 12)

end
