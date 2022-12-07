function foo = generate_random_walk_density(uvt_idx, density_all, N, N_half, density_dir)

    % 3D random walk (steps = 500) to generate a continuous 3D density map

    % uvt_idx:      the index of the random density
    % density_all:  the sum of all the density values, without loss of generality, by default, it is set to 50. The gradient descent step size and regularization parameters are tuned with respect to the case whe density_all is set to 50
    % N:            the odd number N is the size of the sampling Cartesian grid is N by N by N
    % N_half:       half of N, i.e. (N-1)/2
    % density_dir:  a string, the directory where the random densities are saved

    system(['mkdir -p ' density_dir])    % incase the directory was not created yet

    fprintf('random density %d\n', uvt_idx)
    rng(uvt_idx)

    % 3d random walk to determine the mean of the Gaussian source model
    %------------------------------------------------------------------------%
    Num=500;    % the number of random walk steps
    step=1;     % the step size
    Xi=0;       % x starting location
    Yi=0;       % y starting location
    Zi=0;       % z starting location
    X(1:Num)=0;
    Y=X;
    Z=Y;
    X(1)=Xi;
    Y(1)=Yi;
    Z(1)=Zi;
    %------------------------------------------------------------------------%
       
    for i=2:Num
        if(rand(1,1)>=0.5)
            X(i)=X(i-1)+step;
        else
            X(i)=X(i-1)-step;
        end
        if(rand(1,1)<=0.5)
            Y(i)=Y(i-1)+step;
        else
            Y(i)=Y(i-1)-step;
        end
        if(rand(1,1)>=0.5)
            Z(i)=Z(i-1)+step;
        else
            Z(i)=Z(i-1)-step;
        end
             
    end

    s_all = [X' Y' Z'];
    [C ia ic] = unique(s_all,'rows');

    % centeralization
    x_ori = C-mean(s_all);  % shoule be mean(s_all) to account for the weights
    weight_ori = zeros(size(C,1),1);
    for (i=1:size(C,1))
    dist_tmp = sum((s_all-C(i,:)).^2,2);
    weight_ori(i) = length(dist_tmp(dist_tmp==0));
    end
    weight_ori = weight_ori/sum(weight_ori)*density_all; % without loss of generality, normalized the weight vector so that it sums to 50 by default, this stablized the regularization parameters 



    % scale the Gaussian source model so that it fits in the sampling grid
    s_sd_ori = 1;   % the original standard deviation of the Gaussian source model
    N_half_ori = ceil(max(abs(x_ori(:))))+3*s_sd_ori;
    % compute scaling factor to make sure it is within a grid of size 101x101x101
    mut_factor = ceil(0.9*N_half)/N_half_ori; % to be consistent with N=101


    % now scale the density model
    % due to the scaling, there will be a model mismatch between the scaled Gaussian source model and the Cartesian sampling grid
    x_scale = x_ori*mut_factor;   % the scaled Gaussian source mean
    s_sd = s_sd_ori*mut_factor; % the scaled Gaussian source stand deviation


    % compute the final 3d random density sampled on the Cartesian grid
    s_all = x_scale;

    den_3d = repmat(0,[N N N]);
    x_seq = -N_half:N_half;
    y_seq = -N_half:N_half;
    z_seq = -N_half:N_half;   % corresponding to den_3d

    for (s_idx=1:size(s_all,1))
        fprintf('%d\n',s_idx)
        s_tmp = s_all(s_idx,:);
        weight_tmp = weight_ori(s_idx);
        
        
        for (i=(round(s_tmp(1)-3*s_sd)-1):(round(s_tmp(1)+3*s_sd)+1))
            for (j=(round(s_tmp(2)-3*s_sd)-1):(round(s_tmp(2)+3*s_sd)+1))
                for (k=(round(s_tmp(3)-3*s_sd)-1):(round(s_tmp(3)+3*s_sd)+1))
                    den_tmp = weight_tmp * (normcdf(i+0.5,s_tmp(1),s_sd)-normcdf(i-0.5,s_tmp(1),s_sd)) * (normcdf(j+0.5,s_tmp(2),s_sd)-normcdf(j-0.5,s_tmp(2),s_sd)) * (normcdf(k+0.5,s_tmp(3),s_sd)-normcdf(k-0.5,s_tmp(3),s_sd));
                    
                    if (i+N_half+1>=1) && (i+N_half+1<=N)
                    if (j+N_half+1>=1) && (j+N_half+1<=N)
                    if (k+N_half+1>=1) && (k+N_half+1<=N)
                    den_3d(i+N_half+1,j+N_half+1,k+N_half+1) = den_3d(i+N_half+1,j+N_half+1,k+N_half+1) + den_tmp;
                    end
                    end
                    end
                end
            end
        end
    end

    save(strcat(density_dir, '/den_3d_', num2str(uvt_idx), '.mat'), 'den_3d')

end

