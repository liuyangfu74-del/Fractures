function [B, T, S_T_idx, S_B_idx, S_cand_idx] = gmmDoubleThreshold(Nf_ds)
% GMM 三分量定双阈值，T = μ3 + σ3

    Nf_pos = Nf_ds(Nf_ds > 0);
    
    try
        gmm = fitgmdist(Nf_pos, 3, 'Replicates', 5);
        [mu_s, order] = sort(gmm.mu);
        sigma_s = sqrt(squeeze(gmm.Sigma));
        sigma_s = sigma_s(order);
        T = mu_s(3)+2*sigma_s(3);B = mu_s(2)+sigma_s(2);
    end
    
    if T <= B
        T = (B + max(Nf_pos)) / 2;
    end
    
    S_T_idx    = find(Nf_ds > T);
    S_B_idx    = find(Nf_ds > B & Nf_ds <= T);
    S_cand_idx = find(Nf_ds > B);
end