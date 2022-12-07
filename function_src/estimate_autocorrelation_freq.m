function Cl = estimate_partial_correlation(C_FB, info, sampling_points)

load('alpha_ql_maxL_50_maxQ_300.mat')
%load('beta_qs_maxQ_300_lgwt_P_maxS_150.mat')

B_ql = beta_FB_2_bessel_par(300, 150, sampling_points, info.r_cut);

% correct for bandlimit
B_ql = B_ql*sqrt(info.r_cut/0.5);

sz2 = length(sampling_points);

maxQ = length(C_FB)-1;

Cl = compute_autocorrelation_freq_from_alpha_beta(info.maxL, C_FB, alpha_ql, B_ql, maxQ, sz2);
