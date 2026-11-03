function [fracture_mask, S_T_in_cand] = hysteresisTrackingSpatial( ...
    ptCloudDown, Nf_ds, L_i, S_T_idx, S_cand_idx, B)
% 纯空间滞后追踪：只用 visited + Nf > B + dist < L_v 三个判据
%
% 输入：
%   ptCloudDown : 降采样点云
%   Nf_ds       : 降采样点云上的 Nf
%   L_i         : 每个点的自适应边长
%   S_T_idx     : 强点索引（全体点云索引）
%   S_cand_idx  : 候选点索引（全体点云索引）
%   B           : 低阈值
%
% 输出：
%   fracture_mask : M_cand×1 logical，候选集上的裂隙掩码
%   S_T_in_cand   : 强点在候选集中的索引

    % 强点在候选集中的位置
    [~, S_T_in_cand] = ismember(S_T_idx, S_cand_idx);
    
    % 预提取候选集属性
    P_cand  = ptCloudDown.Location(S_cand_idx, :);
    L_cand  = L_i(S_cand_idx);
    Nf_cand = Nf_ds(S_cand_idx);
    M = size(P_cand, 1);
    
    % 初始化
    fracture_mask = false(M, 1);
    fracture_mask(S_T_in_cand) = true;
    
    visited = false(M, 1);
    visited(S_T_in_cand) = true;
    
    merged = false(numel(S_T_in_cand), 1);
    
    kdtree_cand = KDTreeSearcher(P_cand);
    
    % 主循环
    for k = 1:numel(S_T_in_cand)
        if merged(k), continue; end
        
        active = S_T_in_cand(k);
        
        while true
            round_new = [];
            
            for vi = 1:numel(active)
                v = active(vi);
                
                % 在候选集里找 v 的 L_v 邻域
                idx = rangesearch(kdtree_cand, P_cand(v,:), L_cand(v));
                idx = idx{1};
                
                for j = idx
                    if visited(j), continue; end
                    if Nf_cand(j) <= B, continue; end
                    
                    d_vj = norm(P_cand(v,:) - P_cand(j,:));
                    if d_vj >= L_cand(v), continue; end
                    
                    visited(j) = true;
                    fracture_mask(j) = true;
                    round_new(end+1) = j; %#ok<AGROW>
                    
                    % 如果 j 是强点，标记被合并
                    other_k = find(S_T_in_cand == j, 1);
                    if ~isempty(other_k) && other_k ~= k
                        merged(other_k) = true;
                    end
                end
            end
            
            if isempty(round_new), break; end
            active = round_new;
        end
    end
    
    fprintf('合并的强点数：%d / %d\n', sum(merged), numel(S_T_in_cand));
end