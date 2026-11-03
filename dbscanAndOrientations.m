function [cluster_results, idx_crack_clean, labels_clean] = ...
    dbscanAndOrientations(ptCloudDown, idx_crack, L_i, bandwidth_ori)
% DBSCAN 分组 + 逐簇产状 + 段级 D 值过滤
%
% 输入：
%   ptCloudDown   : 降采样点云
%   idx_crack     : 滞后追踪后的裂隙点索引
%   L_i           : 每个点的自适应边长
%   L_max         : 最大裂隙开度（用于 numNeighbors）
%   bandwidth_ori : MeanShift 带宽（度）
%
% 输出：
%   cluster_results : cell，每簇的产状结果
%   idx_crack_clean : 过滤后的裂隙点索引
%   labels_clean    : 每个保留点的簇标签

    %% ---------- DBSCAN 去噪 + 分组 ----------
    fprintf('\n开始 DBSCAN 去噪 ...\n');

    P_crack = ptCloudDown.Location(idx_crack, :);
    kdtree_crack = KDTreeSearcher(P_crack);
    [~, d_nn] = knnsearch(kdtree_crack, P_crack, 'K', 2);
    d_bar_crack = median(d_nn(:, 2));
    eps_db      = 3*d_bar_crack;
    minPts_db   = 24;

%     fprintf('  d_bar = %.4f, eps = %.4f, minPts = %d\n', ...
%             d_bar_crack, eps_db, minPts_db);

    labels_db = dbscan(P_crack, eps_db, minPts_db);
    keep_db   = labels_db > 0;
    idx_crack_clean = idx_crack(keep_db);
    labels_clean    = labels_db(keep_db);

%     fprintf('  DBSCAN 后：%d 点（剔除 %d）\n', ...
%             numel(idx_crack_clean), numel(idx_crack) - numel(idx_crack_clean));

    unique_clusters = unique(labels_clean);
    n_clusters = numel(unique_clusters);
%     fprintf('  共 %d 个裂隙簇\n', n_clusters);

%% ---------- 逐簇产状计算 ----------
fprintf('\n开始逐簇产状计算 ...\n');

cluster_results = cell(n_clusters, 1);

for c = 1:n_clusters
    label_c = unique_clusters(c);
    idx_c   = find(labels_clean == label_c);

    P_c = ptCloudDown.Location(idx_crack_clean(idx_c), :);
    HP_c = idx_crack_clean(idx_c);

    if size(P_c, 1) < 10
%         fprintf('  簇 %d：点数太少（%d），跳过\n', c, size(P_c, 1));
        cluster_results{c} = struct( ...
            'label', label_c, ...
            'n_points', size(P_c, 1), ...
            'azimuths_dips', [], ...
            'DSPCData', []);
        continue;
    end

    % ★ 取出当前簇每个点的 L_i
    L_c = L_i(idx_crack_clean(idx_c));

    % ★ 自适应 numNeighbors：每个点用它自己的 L_i
    kdtree_c = KDTreeSearcher(P_c);
    counts = zeros(size(P_c, 1), 1);
    for i = 1:size(P_c, 1)
        idx_nn = rangesearch(kdtree_c, P_c(i,:), L_c(i));
        counts(i) = numel(idx_nn{1});
    end
    numNeighbors_c = round(median(counts));

    fprintf('  簇 %d：%d 点，numNeighbors = %d\n', ...
            c, size(P_c, 1), numNeighbors_c);

    try
        [azimuths_dips_c, DSPCData_c] = Extractionfractureinformation( ...
            numNeighbors_c, P_c, HP_c, ptCloudDown, bandwidth_ori);

        cluster_results{c} = struct( ...
            'label', label_c, ...
            'n_points', size(P_c, 1), ...
            'azimuths_dips', azimuths_dips_c, ...
            'DSPCData', DSPCData_c);
    catch ME
        fprintf('  簇 %d：产状计算失败（%s）\n', c, ME.message);
        cluster_results{c} = struct( ...
            'label', label_c, ...
            'n_points', size(P_c, 1), ...
            'azimuths_dips', [], ...
            'DSPCData', []);
    end
end

fprintf('产状计算完成。\n');
    %% ---------- 段级 D 值 GMM 过滤 ----------
    fprintf('\n开始段级 D 值过滤 ...\n');

    % 收集所有段的所有 D
    D_all_seg = [];
    seg_to_cluster = [];
    for c = 1:n_clusters
        res = cluster_results{c};
        if isempty(res.DSPCData), continue; end

        n_seg = numel(res.DSPCData);
        for s = 1:n_seg
            D_all_seg = [D_all_seg; res.DSPCData(s).D];
            seg_to_cluster = [seg_to_cluster; c];
        end
    end

    if numel(D_all_seg) < 1
        fprintf('  段数太少（%d），跳过段级滤波。\n', numel(D_all_seg));
    else
        try
            gmm_D = fitgmdist(D_all_seg, 2, 'Replicates', 5);
            labels_seg = cluster(gmm_D, D_all_seg);
            [~, crack_class] = max(gmm_D.mu);
        catch
            warning('GMM 失败，改用 Otsu。');
            D_norm = (D_all_seg - min(D_all_seg)) / ...
                     (max(D_all_seg) - min(D_all_seg) + eps);
            thr = graythresh(D_norm);
            labels_seg = double(D_norm > thr) + 1;
            crack_class = 2;
        end

        % 每簇统计类裂隙段数量
        keep_cluster = false(n_clusters, 1);
        for c = 1:n_clusters
            n_crack_seg = sum(labels_seg(seg_to_cluster == c) == crack_class);
            keep_cluster(c) = (n_crack_seg >= 2);
            fprintf('  簇 %d：类裂隙段 %d 个 → %s\n', ...
                    c, n_crack_seg, ...
                    ternary(keep_cluster(c), '保留', '剔除'));
        end

        % 剔除误判簇
        keep_idx = ismember(labels_clean, unique_clusters(keep_cluster));
        idx_crack_clean = idx_crack_clean(keep_idx);
        labels_clean    = labels_clean(keep_idx);

        fprintf('段级滤波后：%d 簇，%d 点\n', ...
                sum(keep_cluster), numel(idx_crack_clean));
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end