function foo = compute_aligned_model(LN_num, N, N_half, res_cutoff_thd, res_unit, aspire_dir, density_file, input_file, output_file)

    % LN_num:           the largest number of computing threads
    % N:                an odd number, the size of the sampling Cartesian graid is N by N by N
    % N_half:           half of N, i.e. (N-1)/2
    % res_cutoff_thd:       the resolution cutoff threshold based on the FSC curve
    % res_unit:         the length of a single voxel
    % aspire_dir:       a string, the directory where the aspire package is saved
    % density_file:     a string, the location of the density file
    % input_file:       a string, the location of the recovered density file
    % outputfile:       a string, the location of the aligned density file

    
    LN = maxNumCompThreads(LN_num);

    addpath(genpath(aspire_dir))

    den_rec = dlmread(input_file);
    den_3d_rec = repmat(0,[N N N]);
    for (i=1:size(den_rec,1))
        den_3d_tmp = den_rec(i,:);
        den_3d_rec(den_3d_tmp(1)+N_half+1,den_3d_tmp(2)+N_half+1,den_3d_tmp(3)+N_half+1) = den_3d_tmp(4);
    end 

    % compute convolution kernel
    s_sd = 0.8660254;
    s_tmp = [0 0 0];
    cov_kernel = zeros((6*ceil(s_sd)+1),(6*ceil(s_sd)+1),(6*ceil(s_sd)+1));
    for (i=(round(s_tmp(1)-3*ceil(s_sd))-1):(round(s_tmp(1)+3*ceil(s_sd))+1))
    for (j=(round(s_tmp(2)-3*ceil(s_sd))-1):(round(s_tmp(2)+3*ceil(s_sd))+1))
        for (k=(round(s_tmp(3)-3*ceil(s_sd))-1):(round(s_tmp(3)+3*ceil(s_sd))+1))
            if ((i+3*ceil(s_sd)+1 >=1) && (i+3*ceil(s_sd)+1<=(6*ceil(s_sd)+1)))
            if ((j+3*ceil(s_sd)+1 >=1) && (j+3*ceil(s_sd)+1<=(6*ceil(s_sd)+1)))
            if ((k+3*ceil(s_sd)+1 >=1) && (k+3*ceil(s_sd)+1<=(6*ceil(s_sd)+1)))

            den_tmp = normpdf(i,s_tmp(1),s_sd) * normpdf(j,s_tmp(2),s_sd) * normpdf(k,s_tmp(3),s_sd);

            cov_kernel(i+3*ceil(s_sd)+1,j+3*ceil(s_sd)+1,k+3*ceil(s_sd)+1) = cov_kernel(i+3*ceil(s_sd)+1,j+3*ceil(s_sd)+1,k+3*ceil(s_sd)+1) + den_tmp;
            end
            end
            end
        end
    end
    end

    den_3d_rec = convn(den_3d_rec, cov_kernel, 'same');

    load(density_file)

    [Rest,estdx,vol2aligned,reflect]=cryo_align_densities(den_3d,den_3d_rec);

    save(output_file, 'vol2aligned')    % aligned density model

    [resA,fighandle] = plotFSC(den_3d,vol2aligned,res_cutoff_thd,res_unit);

    corr_rec = corr(den_3d(:), vol2aligned(:));
    save(strcat(output_file, '_resA'), 'resA'); % resolution
    save(strcat(output_file, '_corr_rec'), 'corr_rec'); % correlation coefficient


end
