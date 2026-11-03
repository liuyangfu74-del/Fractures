function [B, T, S_T_idx, S_B_idx, S_cand_idx] = maxEntropyThreshold(Nf_ds)
% 四段熵最大，取第二、第三条边界作为 B 和 T

    Nf_pos = Nf_ds(Nf_ds > 0);
    nbins = 128;
    [counts, edges] = histcounts(Nf_pos, nbins);
    centers = (edges(1:end-1) + edges(2:end)) / 2;
    p = counts / sum(counts);
    P = cumsum(p);
    
    best_entropy = -inf;
    best_t = [centers(1), centers(1), centers(end)];
    
    for i = 1:nbins-3
        for j = i+1:nbins-2
            for k = j+1:nbins-1
                p1 = P(i);
                p2 = P(j) - P(i);
                p3 = P(k) - P(j);
                p4 = 1 - P(k);
                if p1 < eps || p2 < eps || p3 < eps || p4 < eps, continue; end
                
                p1_dist = p(1:i) / p1;
                p2_dist = p(i+1:j) / p2;
                p3_dist = p(j+1:k) / p3;
                p4_dist = p(k+1:end) / p4;
                
                H1 = -sum(p1_dist .* log(p1_dist + eps));
                H2 = -sum(p2_dist .* log(p2_dist + eps));
                H3 = -sum(p3_dist .* log(p3_dist + eps));
                H4 = -sum(p4_dist .* log(p4_dist + eps));
                
                H_total = H1 + H2 + H3 + H4;
                
                if H_total > best_entropy
                    best_entropy = H_total;
                    best_t = [centers(i), centers(j), centers(k)];
                end
            end
        end
    end
    
    B = best_t(2);
    T = best_t(3);
    
    S_T_idx    = find(Nf_ds > T);
    S_B_idx    = find(Nf_ds > B & Nf_ds <= T);
    S_cand_idx = find(Nf_ds > B);
end