function foo = generate_mom_init_noisy(density_all, s_sd, num_ref_proj, cumsum_thd_pec, proj_dir, feature_dir, init_dir)

    % density_all:  the sum of all the density values, can be estimated from the spatial radial features
    % s_sd:         the width of Gaussian smooth kernel for the denoised projection images from MFVDM, it is set to 0.5 by default
    % num_ref_proj: the number of reference projections to be used in parallel
    % cumsum_thd_pec:   the hard threshold used to determine candidate locations in the reference projection based on cummulative summation of all the pixels in the reference projection
    % proj_dir:     a string, the directory where the projection images are saved
    % feature_dir:  a string, the directory where the features are saved
    % init_dir:     a string, the directory where the initializations computed from the spatial radial features and a reference image are saved


    load(strcat(proj_dir, '/denoised_images.mat'))
    proj_image_mat = denoised_images;

    N = size(proj_image_mat,1); % N: the odd number N is the size of the sampling Cartesian grid is N by N by N
    N_half = (N-1)/2;   % N_half: half of N, i.e. (N-1)/2

    % smooth and project denoised images by MFVDM
    % maybe smoothing is not needed for low-pass filtered (denoised) images...

    s_tmp = [0 0];
    cov_kernel = zeros((6*ceil(s_sd)+1),(6*ceil(s_sd)+1));
    for (i=(round(s_tmp(1)-3*ceil(s_sd))-1):(round(s_tmp(1)+3*ceil(s_sd))+1))
        for (j=(round(s_tmp(2)-3*ceil(s_sd))-1):(round(s_tmp(2)+3*ceil(s_sd))+1))
            if ((i+3*ceil(s_sd)+1 >=1) && (i+3*ceil(s_sd)+1<=(6*ceil(s_sd)+1)))
                if ((j+3*ceil(s_sd)+1 >=1) && (j+3*ceil(s_sd)+1<=(6*ceil(s_sd)+1)))

                    den_tmp = normpdf(i,s_tmp(1),s_sd) * normpdf(j,s_tmp(2),s_sd);

                    cov_kernel(i+3*ceil(s_sd)+1,j+3*ceil(s_sd)+1) = cov_kernel(i+3*ceil(s_sd)+1,j+3*ceil(s_sd)+1) + den_tmp;
                end
            end
        end
    end

    for (i=1:size(proj_image_mat,3))
        proj_image_mat_tmp = proj_image_mat(:,:,i);
        proj_image_mat_tmp = convn(proj_image_mat_tmp, cov_kernel, 'same');
        proj_image_mat(:,:,i) = proj_image_mat_tmp;
    end

    for (i=1:size(proj_image_mat,3))
        proj_image_mat_tmp = proj_image_mat(:,:,i);
        proj_image_mat_tmp = reshape(ProjectOntoSimplex(proj_image_mat_tmp(:)/density_all), size(proj_image_mat_tmp)) * density_all;
        proj_image_mat(:,:,i) = proj_image_mat_tmp;
    end

    save(strcat(proj_dir, '/denoised_images_smoothed_projected.mat'), 'proj_image_mat')

    % comptue the initializations based on the (uniformly sampled) spatial radial features and a reference projection image

    for (proj_idx=1:num_ref_proj)

        system(['mkdir -p ' init_dir])
        cryo_em_initialization(10, strcat(feature_dir, '/spatial_radial_uniform_feat_noisy'), strcat(proj_dir, '/denoised_images_smoothed_projected.mat'), proj_idx, N, 50, density_all, strcat(init_dir, '/mom_lsq_init_',num2str(proj_idx)))

        load(strcat(init_dir, '/mom_lsq_init_',num2str(proj_idx)))

        % small random perturbation to avoid a symmetric density initialization
        X_init_out = X_init_out.*(1+(rand(size(X_init_out))-0.5)*1);

        den_3d_seq = X_init_out(:);
        den_3d_out = zeros(length(den_3d_seq(den_3d_seq~=0)), 4);
        idx=1;
        for (i=1:N)
            for (j=1:N)
                for (k=1:N)
                    if (X_init_out(i,j,k)~=0)
                        den_3d_out(idx,:) = [i-N_half-1 j-N_half-1 k-N_half-1 X_init_out(i,j,k)];
                        idx = idx+1;
                    end
                end
            end
        end

        den_3d_out(:,4) = den_3d_out(:,4)/sum(den_3d_out(:,4),'all')*density_all;

        %% save the initialization computed from the spatial radial features and the reference projection
        dlmwrite(strcat(init_dir, '/mom_lsq_init_den_3d_',num2str(proj_idx)), den_3d_out, 'delimiter', ' ', 'precision', 12)

        proj_image = proj_image_mat(:,:,proj_idx);
        proj_img_2d = proj_image;

        proj_img_2d_seq = proj_img_2d(:);

        proj_img_2d_seq_norm = sqrt(sum(proj_img_2d_seq.^2));

        proj_img_2d_seq_sort = sort(proj_img_2d_seq, 'descend');
        proj_img_2d_seq_sq_sort = proj_img_2d_seq_sort.^2;
        proj_img_2d_seq_sq_sort_cumsum = cumsum(proj_img_2d_seq_sq_sort);
        cumsum_thd_idx = length(proj_img_2d_seq_sq_sort_cumsum(proj_img_2d_seq_sq_sort_cumsum<=cumsum_thd_pec*proj_img_2d_seq_norm^2));
        cumsum_thd = proj_img_2d_seq_sort(cumsum_thd_idx);
        proj_img_2d(proj_img_2d<cumsum_thd) = 0;

        proj_img_2d_seq = proj_img_2d(:);

        proj_img_2d_out = zeros(length(proj_img_2d_seq(proj_img_2d_seq~=0)),3);
        idx=1;
        for (i=1:N)
            for (j=1:N)
            if (proj_img_2d(i,j)~=0)
                proj_img_2d_out(idx,:) = [i-N_half-1 j-N_half-1 proj_img_2d(i,j)];
                idx = idx+1;
            end
            end
        end

        %% save the reference projection
        dlmwrite(strcat(proj_dir, '/proj_image_2d_out_',num2str(proj_idx)), proj_img_2d_out, 'delimiter', ' ', 'precision', 12);

    end


end
