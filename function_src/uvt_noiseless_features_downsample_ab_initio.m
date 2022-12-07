function foo = uvt_noiseless_features_downsample_ab_initio(LN_num, num_proj_img, N, N_half, downsample_unit, r_cut, l_max, aspire_dir, proj_dir, proj_downsample_dir, feature_downsample_dir)

    % generate noiseless projection images, and estimate the necessary method-of-moments features

    % LN_num:           the largest number of computing threads
    % num_proj_img:     the number of projection images, by default 10000
    % N:                an odd number, the size of the sampling Cartesian graid is N by N by N
    % N_half:           half of N, i.e. (N-1)/2
    % downsample_unit:  the size of the downsampling block, starting from 3 by 3, 5 by 5, 7 by 7, etc. We downsample by summing up all the pixels within the downsampling block
    % r_cut:            frequency cutoff threshold ranging from 0 to 0.5 (or corresponding to 0 to pi), by default it is set to 0.25
    % l_max:            the largest spherical harmonic degree, by default it is set to 10, (2*l_max+1) should be no greater than the number of frequency radius sampling locations due to the rank of autocorrelation feature matrix
    % proj_dir:         a string, the directory where the projection images are saved
    % proj_downsample_dir:  a string, the directory where the downsampled projection images are saved
    % feature_downsample_dir:      a string, the directory where the features are saved


    LN = maxNumCompThreads( LN_num );  % set the largest number of computing threads

    addpath(genpath(aspire_dir))
    initpath

    totTime = tic;

    system(['mkdir -p ' proj_downsample_dir])
    system(['mkdir -p ' feature_downsample_dir])
    load(strcat(proj_dir, '/raw_images.mat'))
    projs = proj_image_mat;
    sigma = 0;

    N_downsample = 2*floor((N_half-(downsample_unit-1)/2)/downsample_unit)+1;
    N_half_downsample = (N_downsample-1)/2;

    downsample_idx_start = (N_half-(downsample_unit-1)/2)-N_half_downsample*downsample_unit + 1; 
    downsample_idx_end = (N_half+1+(downsample_unit-1)/2)+N_half_downsample*downsample_unit;
    downsample_idx_seq = downsample_idx_start:downsample_idx_end;

    projs_downsample = zeros(N_downsample,N_downsample,size(projs,3));
    for (p_idx=1:size(projs,3))
        projs_tmp = projs(:,:,p_idx);
        projs_tmp = projs_tmp(downsample_idx_seq, downsample_idx_seq);
        projs_downsample(:,:,p_idx) = downsample_image(projs_tmp,N_downsample);
    end

    N = N_downsample;
    N_half = N_half_downsample;
    projs = projs_downsample;

    proj_image_mat = projs;
    save(strcat(proj_downsample_dir, '/raw_images.mat'), 'proj_image_mat')

    % frequency sampling along the radial direction upto the frequency cutoff threshold 
    % when r_cut is larger, a denser sampling scheme could be used to improve performance
    k_seq = (0:N_half)/(N) * (r_cut/0.5);   % double sampling
    k_seq(1) = 0.1/(N) * (r_cut/0.5);       % set the first sampling location to a small nonzero value to avoid the singular case
    k_seq_2pi = k_seq * (2*pi); % the c++ program reads k_seq multiplied by 2*pi
    k_idx_seq = 0:(length(k_seq)-1);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% (2*l_max + 1) must be no greater than length(k_seq) due to the rank of the autocorrelation feature matrix
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ((2*l_max+1)>length(k_seq))
        l_max = floor((length(k_seq)-1)/2);
        fprintf('Attention: l_max is reset to %d!!!\n', l_max)

    end

    disp('Done generating projections')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Use steerable PCA to compute autocorrelation %%
    %% the part of code is by Levin et al.          %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    n = num_proj_img;
    downsampleddim = N;

    % **** Change paths! ****

    %%%% CHECK! %%%%
    c = r_cut;
    %R = 50; % Pre calculated % Try volumes of size 330^3, I think the precomputed basis is for R = 150
    R = N_half;

    L0 = downsampleddim;
    n_r = ceil(4*c*R);
    tic_basis = tic;
    [ basis, sample_points ] = precomp_fb( n_r, R, c );
    timing.basis = toc(tic_basis)

    tic_parse = tic;
    def_grp=1;
    ctfs = ones(downsampleddim, downsampleddim, 1);
    ctfs_rad(:,def_grp) = ones(size(sample_points.r));
    timing.parse=toc(tic_parse)

    disp('Done reading')

    %% Prewhitening

    disp('Downsampling')
    tic_whit = tic;

    % Normalize images
    log_message('Normalize background');
    %n = size(projs,1);
    % ********* SKIP BACKGROUND NORMALIZATION ??? *********************

    if sigma == 0
        noise_v_r = 0;
    end
    mean_img = mean(projs, 3);
    %% Now in fourier space
    projs = cfft2(projs);

    num_pool = 20;
    n_im = size(projs,3);

    %% Denoise using ccwf

    % Image is divided by the filter elementwise, so -1
    % This still has some large pixels values
    timing.whit = toc(tic_whit)

    %% Mean estimation and demeaning

    ndef = 1;
    w_f = ones(downsampleddim, downsampleddim);
    w_CTF = ctfs .* repmat(w_f,1,1,ndef);
    regu = 1;
    tic_mean = tic;

    %Solve better conditioned system to get W\mu then get \mu
    ctfid = ones(n, 1);
    mean_image_f = mean_LS( ctfs , ctfid , projs , regu );
    timing.mean = toc(tic_mean)
    mean_image_f = double(mean_image_f);
    projs = double(projs);

    tic_demean = tic;
    [projs] = demean_y_v6( projs , w_CTF , mean_image_f , ctfid );
    projs = real(icfft2(projs));
    timing.demean = toc(tic_demean)

    tic_coeffymu = tic;
    [ coeff_ymu ] = coeff_demean( projs , R , basis , sample_points , num_pool);
    timing.coeffymu = toc(tic_coeffymu)
    [coeff_mean] = coeff_demean( real(icfft2(mean_image_f)) , R , basis , sample_points , 1 );

    %% CTF in new basis: numerical integration

    if mod(L0,2)==1
        w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2)+1 , floor(L0/2)+1:end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
    else
        w_f_rad = interp1( [0:floor(L0/2)] , w_f(floor(L0/2) , floor(L0/2):end) , sample_points.r*((floor(L0/2))/0.5) , 'linear' );
    end

    %% CCWF

    tic_ccwf = tic;
    [ denoised_coeff_ccwf , ~ , ~ , num_eigs, C_FB ] = jobscript_CCWF_cgshrink_jsb( ctfid , w_f_rad , ctfs_rad , basis , sample_points , coeff_mean , coeff_ymu , noise_v_r );
    C_FB{1} = denoised_coeff_ccwf{1}*denoised_coeff_ccwf{1}'/n_im;
    timing.ccwf = toc(tic_ccwf)

    info.r_cut = c;
    info.maxL = l_max; 

    freq_sampling_points = k_seq;   % the frequency sampling locations

    Cl = estimate_autocorrelation_freq(C_FB, info, freq_sampling_points);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End of autocorrelation estimation %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % cor_feat_est_step_2 is the fourier autocorrelation feature
    cor_feat_est_step_2 = zeros([size(Cl{1}) length(Cl)]);
    for (ii=1:length(Cl))
        cor_feat_est_step_2(:,:,ii) = Cl{ii};
    end

    % the autocorrelation computed used Eitan's code needs to be scaled to get the true autocorrelation
    % I understand the first two terms are based on the different formulations between Eitan's paper and steerable PCA, but I don't understand why we need to multiply it by 50, what did they do exactly? this does not seem to be affected by the sum of density values or the sampling grid
    cor_feat_est_step_2 = cor_feat_est_step_2 * ((1/(sqrt(2)*sqrt(pi)) * sqrt(info.r_cut))^2 * 50 *(0.25/r_cut)^2);

    cor_feat_est_step_2_rec = cor_feat_est_step_2;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% compute the spatial radial features %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % note the when l=0, the spherical harmonic coefficients is essentially the first order moments
    info.r_cut = c;
    info.maxL = 0;
    % the frequency sampling points should be determined according to the Gauss-legendre quadrature rule
    [k_seq_lgwt, k_weight_seq_lgwt] = lgwt(N,0,info.r_cut);
    freq_sampling_points = k_seq_lgwt;   % the frequency sampling locations
    k_weight_seq_lgwt_2pi = k_weight_seq_lgwt*(2*pi);

    Cl = estimate_autocorrelation_freq(C_FB, info, freq_sampling_points);
    C0 = Cl{1};
    % the autocorrelation computed used Eitan's code needs to be scaled to get the true autocorrelation
    C0 = C0 * ((1/(sqrt(2)*sqrt(pi)) * sqrt(info.r_cut))^2 * 50 *(0.25/r_cut)^2);

    % compute the first order moments according to cholesky decomposition
    [Ut St] = eig( C0 );
    St_diag = diag(St);
    [St_diag_sort, St_diag_sort_idx] = sort(St_diag, 'descend');
    Ut = Ut(:,St_diag_sort_idx);
    St_diag = St_diag_sort;

    Ut_out = Ut(:,1)*sqrt(St_diag(1));
    Ut_out = real(Ut_out);  % Ut_out might be imaginary

    % note that this should be averaged (or normalized)
    % it should have been divided by (4*pi*k^2), but Y_{00} = sqrt(1/(4*pi)), and A_{00} was not multiplied with k^2 origniall, hence we only need to divide it by sqrt(4*pi)
    radial_feat_est_step_1 = Ut_out/sqrt(4*pi);

    % sampling the radial direction in the spatial domain according to Gauss-lengendre quadrature rule
    % needed to efficiently evaluate the Fourier autocorrelations
    [r_seq_lgwt, r_seq_lgwt_weight]=lgwt(N,0,N_half);
    [r_seq_lgwt_sort, r_seq_lgwt_sort_idx] = sort(r_seq_lgwt,'ascend');
    r_seq_lgwt = r_seq_lgwt_sort;
    r_seq_lgwt_weight = r_seq_lgwt_weight(r_seq_lgwt_sort_idx);

    radial_lgwt_file = [r_seq_lgwt r_seq_lgwt_weight];
    %% the first column contains the Gaussian quadrature sampling nodes
    %% the second column contains the corresponding quadrature weights
    dlmwrite(strcat(feature_downsample_dir, '/spatial_radial_lgwt'), radial_lgwt_file, 'delimiter', ' ', 'precision', 12) 

    r_seq = r_seq_lgwt;
    r_idx_seq = 0:(length(r_seq)-1);

    radial_feat_est_step_2 = zeros(length(r_seq),1);
    k_seq_lgwt_tmp = k_seq_lgwt*(2*pi);
    for (r_idx=1:length(r_seq))
        r = r_seq(r_idx);
        if (r>0)
            k_r_seq_tmp = k_seq_lgwt_tmp*r;
            radial_feat_est_step_2(r_idx) = r*sum(k_seq_lgwt_tmp.*radial_feat_est_step_1.*sin(k_r_seq_tmp).*k_weight_seq_lgwt_2pi);
        else
            k_r_seq_tmp = k_seq_lgwt_tmp*1e-10;
            radial_feat_est_step_2(r_idx) = 1e-10*sum(k_seq_lgwt_tmp.*radial_feat_est_step_1.*sin(k_r_seq_tmp).*k_weight_seq_lgwt_2pi);
        end

    end

    radial_feat_est_out = [r_idx_seq' real(radial_feat_est_step_2)*(2/pi)]; % 2/pi required by formula

    if (sum(r_seq_lgwt_weight.*radial_feat_est_out(:,2))<0)
        radial_feat_est_out(:,2) = -radial_feat_est_out(:,2);
    end

    %% the first column must start from 0 consequtively, it contains the radial indices
    %% the second column contains the corresponding radial integration features, i.e. the first order momoments
    dlmwrite(strcat(feature_downsample_dir, '/spatial_radial_lgwt_feat_noiseless'), radial_feat_est_out, 'delimiter', ' ', 'precision', 12);

    % mapping between the spatial radial indices and the spatial radial value
    radial_map = [r_idx_seq' r_seq];
    %% the first column must start from 0 consequtively, it contains the radial indices
    %% the second column contains the corresonding radial location(values)
    dlmwrite(strcat(feature_downsample_dir, '/spatial_radial_lgwt_map'), radial_map, 'delimiter', ' ', 'precision', 12)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% write fourier correlation features %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cor_feat_est_step_2_rec_sym = cor_feat_est_step_2_rec;
    for (l=0:l_max)
        cor_feat_est_step_2_rec_sym(:,:,l+1) = 0.5*(real(cor_feat_est_step_2_rec(:,:,l+1)) + real(cor_feat_est_step_2_rec(:,:,l+1))');
    end

    %% fourier_ac_decomposition_feat_noiseless_ contains the transformed spherical harmonic coefficients, i.e. A_l \times O_l, where O_l is unknown orthogonal matrices
    %% every line corresponds to the coefficents belonging to a pair of (l, m)-values

    %% l_idx_fourier_ac_decomposition_ contains the corresponding fourier radial indices
    %% every line corresponds to the coefficents belonging to a pair of (l, m)-values

    for (l_out = [ l_max ])

        if (exist(strcat(feature_downsample_dir, '/fourier_ac_decomposition_feat_noiseless_', num2str(l_out))))
            system(['rm ' feature_downsample_dir '/fourier_ac_decomposition_feat_noiseless_' num2str(l_out)]);
        end
        if (exist(strcat(feature_downsample_dir, '/l_idx_fourier_ac_decomposition_', num2str(l_out))))
            system(['rm ' feature_downsample_dir '/l_idx_fourier_ac_decomposition_' num2str(l_out)]);
        end

        for (l=0:l_out)

            [Ut St] = eig(cor_feat_est_step_2_rec_sym(:,:,l+1));
            St_diag = diag(St);
            [St_diag_sort, St_diag_sort_idx] = sort(St_diag, 'descend');    % what if S_diag is negative, it should be all positive, we should not include negative values
            Ut = Ut(:,St_diag_sort_idx);
            St_diag = St_diag_sort;

            %St_diag = St_diag - St_diag(2*l+2);

            l_val = 2*l+1;

            if (l_val>size(Ut,2))
                break;
            end

            Ut_out = Ut(:,1:l_val)*sqrt(diag(St_diag(1:l_val)));
            Ut_out = real(Ut_out);  % Ut_out might be imaginary sometimes

            m_seq = -l:l;
            for (m_idx=1:l_val)
                cor_feat_est_tmp = [l m_seq(m_idx) Ut_out(:,m_idx).'];
                radial_idx_tmp = [l m_seq(m_idx) k_idx_seq];
                dlmwrite(strcat(feature_downsample_dir, '/fourier_ac_decomposition_feat_noiseless_', num2str(l_out)), cor_feat_est_tmp, '-append', 'delimiter', ' ', 'precision', 12)
                dlmwrite(strcat(feature_downsample_dir, '/l_idx_fourier_ac_decomposition_', num2str(l_out)), radial_idx_tmp, '-append', 'delimiter', ' ', 'precision', 12)
            end

        end

    end

    fourier_radial_map = [k_idx_seq' k_seq_2pi'];
    % mapping between the fourier radial indices and the fourier radial value
    %% the first column must start from 0 consequtively, it contains the fourier radial indices
    %% the second column contains the corresponding radial location (values), note that we need to multiply it by 2*pi here, it is uniformly sampled along the Fourier radial direction here
    dlmwrite(strcat(feature_downsample_dir, '/fourier_radial_uniform_map'), fourier_radial_map, 'delimiter', ' ', 'precision', 12)



    % calculate the spatial radial features uniform sampled along the radial direction
    r_seq = (0:N_half);
    r_idx_seq = 0:(length(r_seq)-1);

    radial_feat_est_step_2 = zeros(length(r_seq),1);
    k_seq_lgwt_tmp = k_seq_lgwt*(2*pi);
    for (r_idx=1:length(r_seq))
        r = r_seq(r_idx);
        if (r>0)
            k_r_seq_tmp = k_seq_lgwt_tmp*r;
            radial_feat_est_step_2(r_idx) = r*sum(k_seq_lgwt_tmp.*radial_feat_est_step_1.*sin(k_r_seq_tmp).*k_weight_seq_lgwt_2pi);
        else
            k_r_seq_tmp = k_seq_lgwt_tmp*1e-10;
            radial_feat_est_step_2(r_idx) = 1e-10*sum(k_seq_lgwt_tmp.*radial_feat_est_step_1.*sin(k_r_seq_tmp).*k_weight_seq_lgwt_2pi);
        end

    end

    radial_feat_est_out = [r_idx_seq' real(radial_feat_est_step_2)*(2/pi)]; % 2/pi required by formula

    if (sum(radial_feat_est_out(:,2))<0)
        radial_feat_est_out(:,2) = -radial_feat_est_out(:,2);
    end

    dlmwrite(strcat(feature_downsample_dir, '/spatial_radial_uniform_feat_noiseless'), radial_feat_est_out, 'delimiter', ' ', 'precision', 12);

end
