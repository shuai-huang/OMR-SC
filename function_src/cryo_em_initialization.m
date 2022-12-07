function foo = cryo_em_initialization(LN_num, radial_feat_file, proj_image_file, proj_idx, N, max_cg_ite, density_all, output_file)

    % do spar and exp separately

    sx = N;
    sy = N;
    sz = N;

    LN = maxNumCompThreads( LN_num );  % set the largest number of computing threads

    % generate random sampling pattern accroding to poisson disk
    sx_center = (sx-1)/2+1;
    sy_center = (sy-1)/2+1;
    sz_center = (sz-1)/2+1;

    loc_index = reshape(1:(sx*sy*sz), [sx sy sz]);
    loc_radial_val = zeros(sx*sy*sz,2);
    for (i=1:sx)
        for (j=1:sy)
            for (k=1:sz)
                radial_val_tmp = sqrt((i-sx_center)^2+(j-sy_center)^2+(k-sz_center)^2);
                loc_radial_val(loc_index(i,j,k),1)=loc_index(i,j,k);
                loc_radial_val(loc_index(i,j,k),2)=round(radial_val_tmp);
            end
        end
    end

    sampling_vect = {};
    measurements = [];
    feature_idx = 1;
    % read radial integration feature
    radial_feature_raw = load(radial_feat_file);
    radial_feature = zeros((sx-1)/2+1,2);
    for (i=0:size(radial_feature,1)-1)
        radial_feature(i+1,1) = i;
        radial_feature(i+1,2) = radial_feature_raw(i+1,2);
    end

    for (i=1:size(radial_feature,1))
        sampling_vect{feature_idx} = loc_radial_val(loc_radial_val(:,2)==radial_feature(i,1),1);
        measurements(feature_idx) = radial_feature(i,2);
        feature_idx = feature_idx + 1;
    end

    % read denoised projection image
    load(proj_image_file)
    proj_image = proj_image_mat(:,:,proj_idx);
    proj_image = proj_image/sum(proj_image,'all')*density_all;

    %figure; imshow(proj_image,[])
    for (i=1:sx)
        for (j=1:sy)
            sampling_vect{feature_idx} = squeeze(loc_index(i,j,:));
            measurements(feature_idx) = proj_image(i,j);
            feature_idx = feature_idx+1;
        end
    end
   
    measurements = measurements.';
    feature_num = length(sampling_vect);
    
    % use conjungate gradient descent to find the least square solution with minimum l2 norm
    tol=1e-4;
    mat_sz = [sx sy sz];


    X_init = zeros(sx,sy,sz);
    d_cg = zeros(sx, sy, sz);

    d_cg = d_cg - sx*sy*sz*At_op_3d(A_op_3d(X_init, sampling_vect) - measurements, sampling_vect, mat_sz);

    r_cg = d_cg;
    for (ite=1:max_cg_ite)   % max iteration is 100
        a_cg_n = sum(conj(r_cg).*r_cg, 'all');
        a_cg_d = 0;

        a_cg_d = a_cg_d + sum(sx*sy*sz*(conj(d_cg).*At_op_3d(A_op_3d(d_cg, sampling_vect), sampling_vect, mat_sz)), 'all');

        a_cg_d = real(a_cg_d);
        a_cg = a_cg_n / a_cg_d;
        X_init_pre = X_init;
        X_init = X_init + a_cg * d_cg;

        % enforce positive constraint
        X_init = reshape(ProjectOntoSimplex(X_init(:)/density_all), size(X_init)) * density_all;

        cvg_cg_val = norm(X_init(:)-X_init_pre(:), 'fro')/norm(X_init(:), 'fro');
        fprintf('Ite %d, cvg val %d\n', ite, cvg_cg_val)
        if (cvg_cg_val<tol)
            break;
        end
        r_cg_new = r_cg;
        r_cg_new = r_cg_new - a_cg*sx*sy*sz*At_op_3d(A_op_3d(d_cg, sampling_vect), sampling_vect, mat_sz);
        b_cg = sum(conj(r_cg_new).*r_cg_new, 'all')/sum(conj(r_cg).*r_cg, 'all');
        d_cg = r_cg_new + b_cg*d_cg;
        r_cg = r_cg_new;
    end

    % change it to image domain

    X_init_out = X_init;

    save(output_file, 'X_init_out')


end
