function B_ql = beta_FB_2_bessel_par(maxQ, maxS, freq_sampling_points, r_cut)

% Q is the order of Bessel function
% S is the number of zeros of the Bessel function
% I is the number of zeros of the spherical Bessel function
% r_cut is the frequency radii

maxP = length(freq_sampling_points);

B_ql = zeros((maxQ+1)*maxP*maxS, 1);

% load('FBj_zeros_maxL_100_maxQ_200_maxS_200_maxI_200.mat') % R_FB & R_j
R_FB = zeros(maxQ + 1, maxS);
for q = 0:maxQ
    R_FB(q+1, :) = zerobess('J', q, maxS);
end

parfor jj = 1:(maxQ+1)*maxP*maxS
    % transform index to subscripts:
    [q, p, s] = ind2sub([maxQ+1, maxP, maxS], jj);
    q = q-1;
    
    phi_qs = @(r) besselj(q, r*R_FB(q+1, s)/r_cut)...
        ./(r_cut*sqrt(pi)*abs(besselj(q+1, R_FB(q+1, s)))); % FB func.
    
    B_ql(jj) = phi_qs(freq_sampling_points(p));
end

B_ql = reshape(B_ql, maxQ+1, maxP, maxS);
