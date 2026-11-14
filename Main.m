% ============================================================
% 基于局部裂口几何与方向延伸性的裂隙点云识别
% ============================================================
clear; clc; close all;

% ---------- 1. 读取点云 ----------
ptCloud = pcread('Your_file_path\Test.pcd');

%% ---------- 参数设置 ----------
voxelSize    = 0.05;
Pn           = [-1, -1, 0];
Pn           = Pn / norm(Pn);
q            = 1.5;
alpha        = 8;
thetaThr     = 45;
L_max        = 0.2;

%% ---------- 降采样 ----------
ptCloudDown = pcdownsample(ptCloud, 'gridAverage', voxelSize);

%% ---------- 计算 G_i 和自适应边长 L_i ----------
fprintf('开始计算 G_i 和 L_i ...\n');
[G, L_i] = computeG(ptCloudDown, Pn, q, alpha, L_max);

%% ---------- 对裂隙候选点计算 S_i ----------
fprintf('开始计算 S_i ...\n');
mask = true(size(ptCloudDown.Location, 1), 1);
[S, normals] = computeS(ptCloudDown, L_i, mask, thetaThr, ...
             'Visualize', true, 'nShow', 6, 'rngSeed', 1);

%% ---------- 合成最终响应 N_f ----------
Nf = zeros(size(G));
Nf(mask) = G(mask) .* S(mask);

%% ============================================================
%  分割阶段：四种双阈值方法
% ============================================================

%% ---------- 准备数据 ----------
Nf_ds   = Nf;
N_ds    = numel(Nf_ds);

%% ---------- 四种方法完整跑一遍 ----------
methods = {'GMM', 'Doubleotsu', 'K-means', 'Maximum entropy'};
n_methods = numel(methods);
results = cell(n_methods, 1);
bandwidth_ori = 0.6;   % MeanShift 角度带宽（度）

for m = 1:4
    fprintf('\n========================================\n');
    fprintf('方法 %d：%s\n', m, methods{m});
    fprintf('========================================\n');
    
    %% ① 双阈值
    switch m
        case 1
            [B, T, S_T_idx, S_B_idx, S_cand_idx] = gmmDoubleThreshold(Nf_ds);
        case 2
            [B, T, S_T_idx, S_B_idx, S_cand_idx] = doubleOtsuThreshold(Nf_ds);
        case 3
            [B, T, S_T_idx, S_B_idx, S_cand_idx] = kmeansDoubleThreshold(Nf_ds);
        case 4
            [B, T, S_T_idx, S_B_idx, S_cand_idx] = maxEntropyThreshold(Nf_ds);
    end
    fprintf('  B = %.4f, T = %.4f\n', B, T);
    fprintf('  |S_T| = %d, |S_cand| = %d\n', ...
            numel(S_T_idx), numel(S_cand_idx));
    
    %% ② 滞后追踪
    [fracture_mask, ~] = hysteresisTrackingSpatial( ...
        ptCloudDown, Nf_ds, L_i, S_T_idx, S_cand_idx, B);
    idx_crack = S_cand_idx(fracture_mask);
    fprintf('  滞后追踪：%d 点\n', numel(idx_crack));
    
    %% ③ DBSCAN 分组 + 逐簇产状
    [cluster_results, idx_crack_clean, labels_clean] = ...
        dbscanAndOrientations( ...
            ptCloudDown, idx_crack, L_i, bandwidth_ori);
    
    %% ④ 保存
    results{m} = struct( ...
    'method', methods{m}, ...
    'B', B, 'T', T, ...
    'n_T', numel(S_T_idx), ...
    'n_cand', numel(S_cand_idx), ...
    'n_crack', numel(idx_crack), ...
    'n_clean', numel(idx_crack_clean), ...
    'n_clusters', numel(cluster_results), ...
    'idx_crack_clean', idx_crack_clean, ...
    'labels_clean', labels_clean, ...
    'cluster_results', {cluster_results}, ...
    'S_T_idx', S_T_idx, ...                 % FCP
    'S_SFCP_idx', idx_crack(~ismember(idx_crack, S_T_idx)));  % SFCP
end
fprintf('%-10s %8s %8s %8s %10s %10s %10s\n', ...
        '方法', 'B', 'T', '|S_T|', '追踪后', 'DBSCAN后', '簇数');
for m = 1:4
    r = results{m};
    fprintf('%-10s %8.4f %8.4f %8d %10d %10d %10d\n', ...
            r.method, r.B, r.T, r.n_T, r.n_crack, ...
            r.n_clean, r.n_clusters);
end

%% ---------- Visualization ----------
fontName = 'Times New Roman';
for m = 1:4
    figure;
    hold on;
    
    % Background points (gray)
    bg = setdiff(1:N_ds, results{m}.idx_crack_clean);
    bg = bg(randperm(numel(bg), min(20000, numel(bg))));
    scatter3(ptCloudDown.Location(bg,1), ...
             ptCloudDown.Location(bg,2), ...
             ptCloudDown.Location(bg,3), ...
             1, [0.85 0.85 0.85], 'filled');
    
    % SFCP (blue)
    if ~isempty(results{m}.S_SFCP_idx)
        scatter3(ptCloudDown.Location(results{m}.S_SFCP_idx,1), ...
                 ptCloudDown.Location(results{m}.S_SFCP_idx,2), ...
                 ptCloudDown.Location(results{m}.S_SFCP_idx,3), ...
                 4, [0.0 0.4 0.9], 'filled');
    end
    
    % FCP (orange)
    scatter3(ptCloudDown.Location(results{m}.S_T_idx,1), ...
             ptCloudDown.Location(results{m}.S_T_idx,2), ...
             ptCloudDown.Location(results{m}.S_T_idx,3), ...
             8, [0.95 0.5 0.1], 'filled');
    
    % ===== 最小改动：收紧坐标范围，去掉多余留白 =====
    allX = [ptCloudDown.Location(bg,1); ...
            ptCloudDown.Location(results{m}.S_SFCP_idx,1); ...
            ptCloudDown.Location(results{m}.S_T_idx,1)];
    allY = [ptCloudDown.Location(bg,2); ...
            ptCloudDown.Location(results{m}.S_SFCP_idx,2); ...
            ptCloudDown.Location(results{m}.S_T_idx,2)];
    allZ = [ptCloudDown.Location(bg,3); ...
            ptCloudDown.Location(results{m}.S_SFCP_idx,3); ...
            ptCloudDown.Location(results{m}.S_T_idx,3)];
    
    daspect([1 1 1]);
    margin = 0.05;
    xlim([min(allX)-margin*range(allX), max(allX)+margin*range(allX)]);
    ylim([min(allY)-margin*range(allY), max(allY)+margin*range(allY)]);
    zlim([min(allZ)-margin*range(allZ), max(allZ)+margin*range(allZ)]);
    % ================================================
    
    axis off;
    set(gca, 'LooseInset', get(gca, 'TightInset'));   % 去掉画布留白
    
    title(sprintf('%s\n(FCP = %d, SFCP = %d)', ...
                  results{m}.method, ...
                  numel(results{m}.S_T_idx), ...
                  numel(results{m}.S_SFCP_idx)), ...
          'FontName', fontName, 'FontSize', 16);
    hold off;
    view(285,2);
end
