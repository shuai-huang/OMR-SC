# Orthogonal Matrix Retrieval with Spatial Consensus for 3D Unknown View Tomography

* 3D unknown view tomography (UVT) aims to reconstruct a 3D density map from random 2D projections taken at "unknown" view angles. It arises in applications such as single-particle cryo-electron microscopy (cryo-EM).


* This package contains code files to implement the orthogonal matrix retrieval with spatial consensus (OMR-SC) approach described in the following paper.
```
@article{OMR_SC_UVT,
    author    = {Shuai Huang and Mona Zehni and Ivan Dokmani{\'c} and Zhizhen Zhao},
    title     = {Orthogonal Matrix Retrieval with Spatial Consensus for 3D Unknown-View Tomography},
    journal   = {arXiv preprint},
    volume    = {arXiv:2207.02985},
    year      = {2022},
    url       = {https://arxiv.org/abs/2207.02985},
}
```
If you use this package and find it helpful, please cite the above paper. Thanks :smile:

## Summary
```
    ./function_src      -- This folder contains MATLAB files that are used for feature extraction
    ./omr_sc_src        -- This folder contains C++ files that are used to perform the reconstruction using the proposed OMR-SC approach
    ./library           -- This folder contains the C++ Eigen and zbessel libraries, they are included here for your convenience. Please also install the ASPIRE and MFVDM packages in this folder.
    ./density           -- This folder contains density files saved in MATLAB MAT-file format
    ./proj_images       -- This folder contains random 2D projection images
    ./features          -- This folder contains extracted features
    ./mom_lsq_init      -- This folder contains the initializations computed from the spatial radial features and one reference projection
    ./ab_initio         -- This folder contains the ab initio density models computed using the OMR-SC approach via downsampling
```

## Usage

This package is tested in the UNIX enviroment from Ubuntu, it should also run in other UNIX enviroments like macOS or Debian. 
Detailed comments are within the individual files. You can follow the following steps to run the program.

* Step 0) Download and install the ASPIRE and MFVDM packages:
```
ASPIRE package can be downloaded at https://github.com/PrincetonUniversity/aspire
MFVDM package can be downloaed at https://github.com/frankfyf/MFVDM-cryo

By default, please install both packages to the 'library' folder.
```

* Step 1) Define the directories. Open the file `directory_definitions.m` in the folder `./function_src`, it contains the directory definitions. Go through each definition, and make changes accordingly.

```
density_file = './density/den_3d_1.mat';    % a string, the location of the density file

% update the directories of aspire and MFVDM packages if necessary
aspire_dir = './library/aspire';                    % a string, the directory where the aspire package is saved, updated if necessary
MFVDM_dir = './library/MFVDM-cryo-master';          % a string, the directory where the MFVDM package is saved, update if necessary

proj_noisy_dir = './proj_images/noisy';     % a string, the directory where the projection images are saved
feature_noisy_dir = './features/noisy';     % a string, the directory where the features are saved
init_noisy_dir = './mom_lsq_init/noisy';    % a string, the directory where the initializations computed from the spatial radial features and a reference image are saved

proj_downsample_noisy_dir = './proj_images/noisy_downsample';   % a string, the directory where the downsampled projection images are saved
feature_downsample_noisy_dir = './features/noisy_downsample';   % a string, the directory where the features computed from downsampled projection images are saved
init_downsample_noisy_dir = './mom_lsq_init/noisy_downsample';  % a string, the directory where the initializations computed from the downsampled spatial radial features and a reference image are saved

proj_noiseless_dir = './proj_images/noiseless';     % a string, the directory where the projection images are saved
feature_noiseless_dir = './features/noiseless';     % a string, the directory where the features are saved
init_noiseless_dir = './mom_lsq_init/noiseless';    % a string, the directory where the initializations computed from the spatial radial features and a reference image are saved

proj_downsample_noiseless_dir = './proj_images/noiseless_downsample';   % a string, the directory where the downsampled projection images are saved
feature_downsample_noiseless_dir = './features/noiseless_downsample';   % a string, the directory where the features computed from downsampled projection images are saved
init_downsample_noiseless_dir = './mom_lsq_init/noiseless_downsample';  % a string, the directory where the initializations computed from the downsampled spatial radial features and a reference image are saved

```

* Step 2) Define the variables. Open the file `variable_definitions.m` in the folder `./function_src`, it contains the varialbe definitions. Go through each definition, and make changes accordingly.
```
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
```

* Step 3) Open `MATLAB` and extract the features. Detailed comments are within each individual function file.
```
addpath('./function_src')

% first read the directory and variable definitions
run('./function_src/directory_definitions.m')
run('./function_src/variable_definitions.m')

%%%%%%%%%%%%%%%%%%%%%%%%%
% noisy reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%

% extract noisy features
uvt_noisy_features(LN_num, num_proj_img, N, N_half, r_cut, l_max, snr_level, density_file, aspire_dir, proj_noisy_dir, feature_noisy_dir)

% denoise the noisy projection images
MFVDM_denoise_images_omr_sc(aspire_dir, MFVDM_dir, proj_noisy_dir)

% compute the initializations using the spatial radial features and one reference image, then save the corresponding reference image
generate_mom_init_noisy(density_all, smooth_sd, num_ref_proj, cumsum_thd_pec, proj_noisy_dir, feature_noisy_dir, init_noisy_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%
% the ab initio density model is computed via downsampling
%%%%%%%%%%%%%%%%%%%%%%%%%

% extract noisy features from downsampled projection images
uvt_noisy_features_downsample_ab_initio(LN_num, num_proj_img, N, N_half, downsample_unit, r_cut, l_max_downsample, aspire_dir, proj_noisy_dir, proj_downsample_noisy_dir, feature_downsample_noisy_dir)

% denoise the downsampled noisy projection images
MFVDM_denoise_images_omr_sc( aspire_dir, MFVDM_dir, proj_downsample_noisy_dir)

% compute the initializations using the downsampled spatial radial features and one downsampled reference image, then save the corresponding downsampled reference image
generate_mom_init_noisy(density_all, smooth_sd, num_ref_proj, cumsum_thd_pec, proj_downsample_noisy_dir, feature_downsample_noisy_dir, init_downsample_noisy_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%
% noiseless reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%

% extract noiseless feature
uvt_noiseless_features(LN_num, num_proj_img, N, N_half, r_cut, l_max, snr_level, density_file, aspire_dir, proj_noiseless_dir, feature_noiseless_dir)

% compute the initializations using the spatial radial features and one reference image, then save the corresponding reference image
generate_mom_init_noiseless(density_all, num_ref_proj, cumsum_thd_pec, proj_noiseless_dir, feature_noiseless_dir, init_noiseless_dir)

%%%%%%%%%%%%%%%%%%%%%%%%%
% the ab initio density model is computed via downsampling
%%%%%%%%%%%%%%%%%%%%%%%%%

% extract noiseless features from downsampled projection images
uvt_noiseless_features_downsample_ab_initio(LN_num, num_proj_img, N, N_half, downsample_unit, r_cut, l_max_downsample, aspire_dir, proj_noiseless_dir, proj_downsample_noiseless_dir, feature_downsample_noiseless_dir)

% compute the initializations using the downsampled spatial radial features and one downsampled reference image, then save the corresponding downsampled reference image
generate_mom_init_noiseless(density_all, num_ref_proj, cumsum_thd_pec, proj_downsample_noiseless_dir, feature_downsample_noiseless_dir, init_downsample_noiseless_dir)

```

* Step 4) Exit `MATLAB`. Then go to the `omr_sc_src` folder, and compile the C++ program
```
g++ -std=c++17 -fopenmp -O3 -o main ./main.cpp ./DataReader.cpp ./UVT.cpp ./PGD.cpp ./global.cpp -I ../library/eigen-3.3.9 -I ../library/zbessel-master -pthread
```

* Step 5) The `main` program takes 11 arguments, it can be run through a script like `run_main_noisy.sh`
```
radial_file                --- the spatial radial features at the spatial radial sampling locations determined according to the "Legendre-Gauss Quadrature Weights and Nodes"
radial_map_file            --- the spatial radial sampling locations determined according to the "Legendre-Gauss Quadrature Weights and Nodes"
radial_lgwt_weight_file    --- the "Legendre-Gauss Quadrature Weights and Nodes"
proj_file                  --- the reference projection image, i.e. the linear projection features
fourier_cor_file           --- the rotated spherical harmonic coefficients obtained by decompositon the Fourier autocorrelation feature matrix through cholesky decomposition
fourier_radial_map_file    --- the Fourier radial sampling locations
fourier_l_idx_file         --- the index file used to map the rotated spherical harmonic coefficients to the (l,m,k) triplet
option_file                --- the option file;
init_file                  --- the initialization file
output_file                --- the recovered density file
obj_file                   --- the objection function file

```


* Step 6) Set the base options files (that end with the suffix "_base") in the four folders `options_downsample_noisy`, `options_noisy`, `options_downsample_noiseless`, and `options_noiseless`. Take the base option file "options_ptn_10_base" in `options_noisy` for example, it is the option file when the largest spherical harmonic degree "l_max" is 10. Each option occupies on line, it starts with the option name, then a space, followed by a numeric value. Please adjust each option accordingly. Note that the comments within the paratheses are added here to explain each option and should not be included in the actual option file.
```

num_thread 10               (the number of computing threads)
r_max 50                   （the maximum radius）
L_max 10                   （the maximum spherical harmonic degree）
num_fourier_l_idx_seq 121   (the number of spherical harmonic (l,m) pairs given L_max)
num_fourier_cor_seq 121     (the number of spherical harmonic (l,m) pairs given L_max)
N 101                       (an odd number, the size of the sampling Cartesian graid is N by N by N)
num_radial 101              (the number of spatial radial sampling locations)
num_fourier_radial 51       (the number of fourier radial sampling locations)
init_type 2                 (the initialization type, "1" for the random initialization; "2" for the prespecified initialization)
max_alt_ite 1000            (the maximum alternating iterations)
max_ite 10                  (the maximum project gradient descent iterations in each inner loop)
density_all 50              (the sum of all the density values, can be estimated from the spatial radial features accurately)
voxel_std 0.8660254         (the standard deviation in the bump function which is chosen to be the isotropic Gaussian function, it is set to sqrt{3}/2 empirically)
bkt_rate 0.95               (the backtracking rate of the gradient descent step)
step_ori 1e-6               (the initial gradient descent step)
step_thd 1e-8               (the minimum gradient descent step)
step_max 1e5                (the maximum gradient descent step)
cvg_thd 1e-6                (the convergence threshold)
```
There are three more options that are set by the script `run_main_noisy.sh`
```
radial_weight TBD           (the regularization parameter for the spatial radial features)
proj_weight TBD             (the regularization parameter for the linear projection features, i.e. the reference image)
num_proj  TBD               (the number of linear projection features, it varies for each reference projection)
```

* Step 7) Take the script `run_main_noisy.sh` for example, it takes four arguments
```
l   --- the maximum spherical harmonic degree, i.e. the previously defined l_max
p   --- the regularization parameter for the spatial radial features and the linear projection features
r   --- the downsample unit
m   --- the index for the reference projection, here we use 10 reference projections for parallel reconstructions, hence m varies from 1 to 10
```
The 4 directories `feature_dir`, `proj_dir`, `option_dir`, `init_dir`, `results_dir` need to be set properly as follows.
```
feature_dir="../features/noisy"
proj_dir="../proj_images/noisy"
option_dir="./options_noisy"
init_dir="../ab_initio/noisy"
results_dir="./results_noisy"
```

* Step 8) Recover the ab initio density models using the features extracted from downsampled projections. If you have access to a computer cluster, you can modify the following scripts to submit the jobs to the cluster using the `qsub` command, so that the jobs can be run in parallel.
```
./run_job_main_noisy_downsample.sh
./run_job_main_noiseless_downsample.sh
```

* Step 9) Go to the parent directory, open `MATLAB`, and upsample the ab initio models.
```
addpath('./function_src')

% first read the directory and variable definitions
run('./function_src/directory_definitions.m')
run('./function_src/variable_definitions.m')

% define the ab initio directories
ab_initio_noisy_dir = './ab_initio/noisy';
ab_initio_noiseless_dir = './ab_initio/noiseless';
low_res_noisy_dir = './omr_sc_src/results_noisy_downsample';
low_res_noiseless_dir = './omr_sc_src/results_noiseless_downsample';

% noisy ab initio modeling
% idx is the index for the reference projection
for (idx = 1:num_ref_proj)
    low_res_file = strcat(low_res_noisy_dir, '/output_', num2str(l_max_downsample), '_', num2str(reg_par), '_', num2str(downsample_unit), '_',num2str(idx));
    ab_initio_file = strcat(ab_initio_noisy_dir, '/init_', num2str(l_max_downsample), '_', num2str(reg_par), '_', num2str(downsample_unit), '_',num2str(idx));
    
    uvt_convert_ab_initio(density_all, N, N_half, downsample_unit, low_res_file, ab_initio_file)
end

% noiseless ab initio modeling
% idx is the index for the reference projection
for (idx = 1:num_ref_proj)
    low_res_file = strcat(low_res_noiseless_dir, '/output_', num2str(l_max_downsample), '_', num2str(reg_par), '_', num2str(downsample_unit), '_',num2str(idx));
    ab_initio_file = strcat(ab_initio_noiseless_dir, '/init_', num2str(l_max_downsample), '_', num2str(reg_par), '_', num2str(downsample_unit), '_',num2str(idx));

    uvt_convert_ab_initio(density_all, N, N_half, downsample_unit, low_res_file, ab_initio_file)
end
```

* Step 10) Exit `MATLAB`. Go the folder `omr_sc_src`, and reconstruct the final density maps. If you have access to a computer cluster, you can modify the following scripts to submit the jobs to the cluster using the `qsub` command, so that the jobs can be run in parallel.
```
./run_job_main_noisy.sh
./run_job_main_noiseless.sh
```

* Step 11) Go to the parent directory: `cd ..`, open `MATLAB`, and align the recovered densities with the groundtruth density.
```
% Among the parallel reconstructions, select the one that minimizes the MSE of the autocorrelation features

addpath('./function_src')

% first read the directory and variable definitions
run('./function_src/directory_definitions.m')
run('./function_src/variable_definitions.m')

% noisy reconstruction
obj = zeros(num_ref_proj,1);
for (idx = 1:num_ref_proj)
	obj_tmp = dlmread( strcat('./omr_sc_src/results_noisy/obj_', num2str(l_max), '_', num2str(reg_par), '_', num2str(downsample_unit), '_', num2str(idx)) );
	obj(idx) = obj_tmp(1);
end

[obj_min, idx_min] = min(obj);

idx = idx_min;
density_rec_file = strcat('./omr_sc_src/results_noisy/output_', num2str(l_max), '_', num2str(reg_par), '_', num2str(downsample_unit), '_', num2str(idx));
density_rec_aligned_file = strcat('./omr_sc_src/results_noisy/output_', num2str(l_max), '_', num2str(reg_par), '_', num2str(downsample_unit), '_', num2str(idx), '_aligned');
compute_aligned_model(LN_num, N, N_half, res_cutoff_thd, res_unit, aspire_dir, density_file, density_rec_file, density_rec_aligned_file)

load(strcat(density_rec_aligned_file, '_corr_rec.mat'))
fprintf('Autocorrelation coefficient: %d\n', corr_rec)


% noiseless reconstruction
obj = zeros(num_ref_proj,1);
for (idx = 1:num_ref_proj)
	obj_tmp = dlmread( strcat('./omr_sc_src/results_noiseless/obj_', num2str(l_max), '_', num2str(reg_par), '_', num2str(downsample_unit), '_', num2str(idx)) );
	obj(idx) = obj_tmp(1);
end

[obj_min, idx_min] = min(obj);

idx = idx_min;
density_rec_file = strcat('./omr_sc_src/results_noiseless/output_', num2str(l_max), '_', num2str(reg_par), '_', num2str(downsample_unit), '_', num2str(idx));
density_rec_aligned_file = strcat('./omr_sc_src/results_noiseless/output_', num2str(l_max), '_', num2str(reg_par), '_', num2str(downsample_unit), '_', num2str(idx), '_aligned');
compute_aligned_model(LN_num, N, N_half, res_cutoff_thd, res_unit, aspire_dir, density_file, density_rec_file, density_rec_aligned_file)

load(strcat(density_rec_aligned_file, '_corr_rec.mat'))
fprintf('Autocorrelation coefficient: %d\n', corr_rec)

```
