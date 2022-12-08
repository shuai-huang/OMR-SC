LN_num = 10;            % the largest number of computing threads
num_proj_img = 10000;   %  the number of projection images, by default 10000
N = 101;                % an odd number, the size of the sampling Cartesian graid is N by N by N
N_half = 50;            % half of N, i.e. (N-1)/2
r_cut = 0.25;           % frequency cutoff threshold ranging from 0 to 0.5 (or corresponding to 0 to pi), by default it is set to 0.25
l_max = 10;             % the largest spherical harmonic degree, by default it is set to 10, (2*l_max+1) should be no greater than the rank of the autocorrelation feature matrix
snr_level = 0.1;        % the signal to noise ratio, i.e. the power of signal divided by the power of noise, by default it is set to 0.1

num_ref_proj = 10;      % the number of reference projections to be used in parallel
density_all = 50;       % the sum of all the density values, can be estimated from the spatial radial features. Without loss of generality, all the densities are normalized so that density_all is 50 here, so that the gradient descent step size and the regularization parameter are relatively stable
smooth_sd = 0.5;        % the smoothing kernel width for denoised projection images (relative to a pixel)
cumsum_thd_pec = 0.99;  % the hard threshold used to determine candidate locations in the reference projection based on the normalized mean squared error between the threshoded reference projection and the original reference projection, "0.99" means the normalized mean squared error would be "0.01". 

downsample_unit = 3;    % the size of the downsampling block, starting from 3 by 3, 5 by 5, 7 by 7, etc. We downsample by summing up all the pixels within the downsampling block
l_max_downsample = 8;   % the largest spherical harmonic degree for downsampled density model, (2*l_max+1) should be no greater than the rank of the autocorrelation feature matrix

res_cutoff_thd = 0.5;   % the resolution cutoff threshold based on the FSC curve
res_unit = 1;           % the length of a single voxel, needs to be updated for different density models

reg_par = 100;          % the regularization weights for the spatial radial features and the linear projection features (the reference image), sometimes setting the regularization weight for spatial radial features to 0 might help, othertimes it might not though... the decomposition of autocorrelation features with l=0 already enforces spatial radial features
