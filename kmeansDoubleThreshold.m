function [B, T, S_T_idx, S_B_idx, S_cand_idx] = kmeansDoubleThreshold(Nf_ds)
% K=4 聚类，取第二、第三条类间边界作为 B 和 T

    Nf_pos = Nf_ds(Nf_ds > 0);
    [labels, centers] = kmeans(Nf_pos, 4, 'Replicates', 5);
    [~, order] = sort(centers);
    
    x1 = Nf_pos(labels == order(1));
    x2 = Nf_pos(labels == order(2));
    x3 = Nf_pos(labels == order(3));
    x4 = Nf_pos(labels == order(4));
    
    B = (max(x1) + min(x2)) / 2;   % 第一条边界 → 不是 B
    T = (max(x2) + min(x3)) / 2;   % 第二条边界
    
    % 等等，B 应该取第几条边界？
    % 如果四类对应：背景、过渡低、过渡高、裂隙
    % B 应该取类 1 和类 2 的边界
    % T 应该取类 3 和类 4 的边界
    
    B = (max(x1) + min(x2))/2;
    T = (max(x3) + min(x4))/2;
    
    S_T_idx    = find(Nf_ds > T);
    S_B_idx    = find(Nf_ds > B & Nf_ds <= T);
    S_cand_idx = find(Nf_ds > B);
end