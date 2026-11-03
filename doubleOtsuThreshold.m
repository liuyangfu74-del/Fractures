function [B, T, S_T_idx, S_B_idx, S_cand_idx] = doubleOtsuThreshold(Nf_ds)
% 四分类 Otsu，取第二、第三条类间边界作为 B 和 T

    Nf_pos = Nf_ds(Nf_ds > 0);
    nbins = 256;
    [counts, edges] = histcounts(Nf_pos, nbins);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    p = counts / sum(counts);
    
    omega = cumsum(p);
    mu = cumsum(centers .* p);
    mu_total = mu(end);
    
    best_sigma = -inf;
    best_t = [centers(1), centers(1), centers(end)];  % 三个阈值
    
    for i = 1:nbins-3
        for j = i+1:nbins-2
            for k = j+1:nbins-1
                w1 = omega(i);
                w2 = omega(j) - omega(i);
                w3 = omega(k) - omega(j);
                w4 = 1 - omega(k);
                if w1 < eps || w2 < eps || w3 < eps || w4 < eps, continue; end
                
                m1 = mu(i) / w1;
                m2 = (mu(j) - mu(i)) / w2;
                m3 = (mu(k) - mu(j)) / w3;
                m4 = (mu_total - mu(k)) / w4;
                m_all = w1*m1 + w2*m2 + w3*m3 + w4*m4;
                
                sigma_b = w1*(m1-m_all)^2 + w2*(m2-m_all)^2 + ...
                          w3*(m3-m_all)^2 + w4*(m4-m_all)^2;
                if sigma_b > best_sigma
                    best_sigma = sigma_b;
                    best_t = [centers(i), centers(j), centers(k)];
                end
            end
        end
    end
    
    B = best_t(2);   % 第二条边界
    T = best_t(3);   % 第三条边界
    
    S_T_idx    = find(Nf_ds > T);
    S_B_idx    = find(Nf_ds > B & Nf_ds <= T);
    S_cand_idx = find(Nf_ds > B);
end