function [G, L_i] = computeG(ptCloudDown, Pn, q, alpha, L_max)
% 计算每个点的局部裂口几何响应 G_i 和自适应边长 L_i
%
% 自适应边长 L_i：取 C_i(L) + D_i(L) 峰值对应的尺度
% 凹陷响应 D_i：切平面内质心到几何中心的距离 / L
% 最终 G_i：在 L_i 下算 C、D，全体点统一稳健归一化后合成

    N = size(ptCloudDown.Location, 1);
    G = zeros(N, 1);
    L_i = zeros(N, 1);

    % ---------- 平均点间距 ----------
    kdtree = KDTreeSearcher(ptCloudDown.Location);
    [~, dist] = knnsearch(kdtree, ptCloudDown.Location, 'K', 2);
    d_bar = median(dist(:, 2));
    L_min = 2 * d_bar;
     
    % ---------- 边长序列 ----------
    L_seq = L_min;
    while L_seq(end) * q <= L_max
        L_seq(end+1) = L_seq(end) * q;
    end
    if L_seq(end) < L_max
        L_seq(end+1) = L_max;
    end
    K = length(L_seq);

    fprintf('  平均点间距 d_bar = %.4f m\n', d_bar);
    fprintf('  边长序列 L_min ~ L_max: %.4f ~ %.4f m，共 %d 档\n', L_min, L_max, K);

        % ---------- 第一遍：算每档边长下的 C_i、D_i ----------
    C_all   = nan(N, K);
    D_all   = nan(N, K);

    for i = 1:N
        p_i = ptCloudDown.Location(i, :);
        for k = 1:K
            idx = rangesearch(kdtree, p_i, L_seq(k));
            idx = idx{1};
            if length(idx) < 4
                continue;
            end
            pts = ptCloudDown.Location(idx, :);
            c_i = mean(pts, 1);

            % PCA（不去质心）
            [coeff, ~, latent] = pca(pts);
            C_all(i, k) = min(latent) / (sum(latent) + eps);

            % 切平面基向量
            e1 = coeff(:, 1);
            e2 = coeff(:, 2);

            % 只投影 p_i 和 c_i
            p_proj1 = dot(p_i - p_i, e1);    % = 0
            p_proj2 = dot(p_i - p_i, e2);    % = 0

            c_proj1 = dot(c_i - p_i, e1);
            c_proj2 = dot(c_i - p_i, e2);

            % 投影质心距
            d_proj = sqrt((c_proj1 - p_proj1)^2 + (c_proj2 - p_proj2)^2);

            % ---------- 三维质心距 ----------
            d_3d = norm(c_i - p_i);

            % ---------- 凹陷响应：法向分量（勾股）----------
            D_all(i, k) = sqrt(max(d_3d^2 - d_proj^2, 0));
        end
    end

    % ============================================================
    %  阶段 1：用 C + D 峰值找 L_i
    % ============================================================
    n_invalid = 0;

    for i = 1:N
        P_curve = nan(K, 1);
        for k = 1:K
            if isnan(C_all(i,k)) || isnan(D_all(i,k))
                continue;
            end
            P_curve(k) = C_all(i,k) + D_all(i,k);
        end

        valid = ~isnan(P_curve);
        if sum(valid) < 2
            L_i(i) = L_max;
            n_invalid = n_invalid + 1;
            continue;
        end

        [~, idx_max] = max(P_curve(valid));
        L_valid = L_seq(valid);
        L_i(i) = L_valid(idx_max);
    end

    % ============================================================
    %  阶段 2：在 L_i 下算最终 G_i（全体点统一归一化）
    % ============================================================
    C_final = nan(N, 1);
    D_final = nan(N, 1);
    Dot_final = nan(N, 1);

    for i = 1:N
        L = L_i(i);
        if L <= 0
            continue;
        end
        idx = rangesearch(kdtree, ptCloudDown.Location(i,:), L);
        idx = idx{1};
        if length(idx) < 10
            continue;
        end
        pts = ptCloudDown.Location(idx, :);
        c_i = mean(pts, 1);

        [coeff, ~, latent] = pca(pts - c_i);
        C_final(i) = min(latent) / (sum(latent) + eps);

        % ---------- 投影到切平面 ----------
        e1 = coeff(:, 1);
        e2 = coeff(:, 2);
        rel = pts - ptCloudDown.Location(i,:);
        proj1 = rel * e1;
        proj2 = rel * e2;

        % 几何中心
        m1 = median(proj1);
        m2 = median(proj2);

        % 质心投影
        c_proj1 = mean(proj1);
        c_proj2 = mean(proj2);

        % 凹陷响应
        d_center = sqrt((c_proj1 - m1)^2 + (c_proj2 - m2)^2);
        D_final(i) = d_center / L;

        % 方向约束（用三维质心向量）
        Dot_final(i) = dot(c_i - ptCloudDown.Location(i,:), Pn);
    end

    % ---------- 全体点统一稳健归一化 + 截断 ----------
    C_final_tilde = (C_final - median(C_final, 'omitnan')) / ...
                    max(iqr(C_final(~isnan(C_final))), eps);
    D_final_tilde = (D_final - median(D_final, 'omitnan')) / ...
                    max(iqr(D_final(~isnan(D_final))), eps);
    C_final_tilde = max(min(C_final_tilde, 5), -5);
    D_final_tilde = max(min(D_final_tilde, 5), -5);

    % ---------- 最终 G_i（加法 + 方向约束）----------
    sig_final = 1 ./ (1 + exp(-alpha * Dot_final ./ L_i));
    G = (C_final_tilde) .* sig_final;
    G = max(G, 0);
    G(isnan(G)) = 0;

end

function plotLidistribution(ptCloudDown, L_i, L_seq, C_all, D_all)
% Visualize the adaptive side length L_i.

    fontName = 'Times New Roman';
    fs_label = 18;
    fs_tick  = 16;
    fs_title = 18;

        % 选一个代表点：L_i 最接近中位数
    L_med = median(L_i(L_i >0));
    [~, idx_show] = min(abs(L_i - L_med));

    %% Response curve
    P_curve  = C_all(idx_show, :) + D_all(idx_show, :);
    figure('Name', 'Response curve', ...
           'Units', 'centimeters', ...
           'Position', [10 10 18 6]);
    plot(L_seq, P_curve, 'b-', 'LineWidth', 2);
    xlabel('Side length L (m)', 'FontName', fontName, 'FontSize', fs_label);
    ylabel('C_i(L) + D_i(L)', 'FontName', fontName, 'FontSize', fs_label);
    set(gca, 'FontName', fontName, 'FontSize', fs_tick);
    title(sprintf('Local geometric response vs side length'), ...
          'FontName', fontName, 'FontSize', fs_title);
    grid off;
    box on;
end