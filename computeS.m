function [S, normals] = computeS(ptCloudDown, L_i, mask, thetaThr, varargin)
% 计算每个裂隙候选点的方向延伸性响应 S_i 和法向量 normals
%
% 输入：
%   ptCloudDown : 降采样点云
%   L_i         : 每个点的自适应边长（N×1）
%   mask        : 裂隙候选掩码（N×1 logical）
%   thetaThr    : 方向一致性阈值（度），本文取 45
%   varargin    : 可选参数
%       'Visualize', true/false  （默认 false）
%       'nShow', 数值            （默认 6）
%       'rngSeed', 数值          （默认 1）
%
% 输出：
%   S         : 方向延伸性响应，二值 {0, 1}
%   normals   : N×3，每个点的局部法向量（来自中心格子 PCA 的最小特征值方向）
%               非候选点或无有效邻域的点，法向量为 [0 0 1]

    %% ---------- 解析可选参数 ----------
    p = inputParser;
    addParameter(p, 'Visualize', false, @islogical);
    addParameter(p, 'nShow', 6, @isnumeric);
    addParameter(p, 'rngSeed', 1, @isnumeric);
    parse(p, varargin{:});
    doVis   = p.Results.Visualize;
    nShow   = p.Results.nShow;
    rngSeed = p.Results.rngSeed;

    %% ---------- 初始化 ----------
    N = size(ptCloudDown.Location, 1);
    S = zeros(N, 1);
    normals = zeros(N, 3);
    normals(:, 3) = 1;                    % 默认法向 (0,0,1)
    kdtree = KDTreeSearcher(ptCloudDown.Location);

    n_center_empty = 0;

    %% ---------- 主循环 ----------
    for i = 1:N
        if ~mask(i) || L_i(i) <= 0
            continue;
        end

        p_i = ptCloudDown.Location(i, :);
        L = L_i(i);

        % 邻域（覆盖九宫格外接圆）
        idx = rangesearch(kdtree, p_i, 1.5 * L * sqrt(2));
        idx = idx{1};
        if length(idx) < 10
            continue;
        end
        pts = ptCloudDown.Location(idx, :);
        rel = pts - p_i;

        % 整体邻域 PCA 定切平面
        [coeff, ~, ~] = pca(rel);
        e1 = coeff(:, 1);
        e2 = coeff(:, 2);

        % 投影到切平面
        proj1 = rel * e1;
        proj2 = rel * e2;

        % 中心格子
        center_in = find(abs(proj1) < L/2 & abs(proj2) < L/2);
        if length(center_in) < 10
            n_center_empty = n_center_empty + 1;
            continue;
        end

        % 中心格子 PCA
        [coeff_c, ~, latent_c] = pca(rel(center_in, :));
        d0 = coeff_c(:, 1);

        % ★ 输出法向量：中心格子的最小特征值方向
        normals(i, :) = coeff_c(:, 2)';
        % 归一化（保险）
        nn = norm(normals(i, :));
        if nn > eps
            normals(i, :) = normals(i, :) / nn;
        else
            normals(i, :) = [0 0 1];
        end

        % 遍历周围 8 个格子
        offsets = [-L, 0, L];
        n_coh = 0;

        for a = 1:3
            for b = 1:3
                if a == 2 && b == 2
                    continue;
                end
                c1 = offsets(a);
                c2 = offsets(b);
                in_cell = find(proj1 > c1 - L/2 & proj1 < c1 + L/2 & ...
                               proj2 > c2 - L/2 & proj2 < c2 + L/2);
                if length(in_cell) < 10
                    continue;
                end
                [coeff_j, ~, ~] = pca(rel(in_cell, :));
                d_j = coeff_j(:, 1);

                cos_theta = abs(dot(d_j, d0));
                cos_theta = min(max(cos_theta, -1), 1);
                theta_j = acosd(cos_theta);
                if theta_j < thetaThr
                    n_coh = n_coh + 1;
                end
            end
        end

        if n_coh >= 1
            S(i) = 1;
        else
            S(i) = 0;
        end
    end

    %% ---------- 输出统计 ----------
    fprintf('  S_i 计算完成。S_i = 1 的点数：%d / %d\n', sum(S), sum(mask));
    fprintf('  中心格子点数不足而跳过的点数：%d\n', n_center_empty);

%     %% ---------- 可视化 ----------
%     if doVis
%         visualizeS(ptCloudDown, L_i, mask, S, thetaThr, nShow, rngSeed);
%     end
end
%% ============================================================
%  可视化函数
% ============================================================
function visualizeS(ptCloudDown, L_i, mask, S, thetaThr, nShow, rngSeed)

    N = length(S);

    %% ---------- 1. S_i 二值分布 ----------
    figure('Name', 'S_i 二值分布');
    histogram(S(mask), [-0.5 0.5 1.5]);
    xticks([0 1]);
    xticklabels({'S_i=0', 'S_i=1'});
    xlabel('S_i'); ylabel('频数');
    title('S_i 二值分布（仅裂隙候选点）');

    %% ---------- 2. S_i 三维分布 ----------
    figure('Name', 'S_i 三维分布');
    scatter3(ptCloudDown.Location(mask,1), ...
             ptCloudDown.Location(mask,2), ...
             ptCloudDown.Location(mask,3), 5, S(mask), 'filled');
    axis equal;
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title('S_i 三维分布（仅裂隙候选点）');
    colorbar; colormap(jet);
    caxis([0, 1]);
    view(45, 30);

    %% ---------- 4. 统计量 ----------
    fprintf('  S_i = 1 占比（裂隙候选内）：%.2f%%\n', ...
            100*sum(S(mask))/sum(mask));

    %% ---------- 5. 单点九宫格可视化 ----------
    rng(rngSeed);
    candidates = find(S == 1 & mask);
    if isempty(candidates)
        fprintf('  没有 S_i = 1 的点，跳过单点可视化。\n');
        return;
    end
    nShow = min(nShow, length(candidates));
    idx_show = candidates(randperm(length(candidates), nShow));

    kdtree = KDTreeSearcher(ptCloudDown.Location);

    for t = 1:nShow
        i = idx_show(t);
        p_i = ptCloudDown.Location(i, :);
        L = L_i(i);

        idx = rangesearch(kdtree, p_i, 1.5 * L * sqrt(2));
        idx = idx{1};
        pts = ptCloudDown.Location(idx, :);
        rel = pts - p_i;

        [coeff, ~, ~] = pca(rel);
        e1 = coeff(:, 1);
        e2 = coeff(:, 2);
        n_local = coeff(:, 3);

        proj1 = rel * e1;
        proj2 = rel * e2;

        center_in = find(abs(proj1) < L/2 & abs(proj2) < L/2);
        if length(center_in) < 3
            continue;
        end
        [coeff_c, ~, ~] = pca(rel(center_in, :));
        d0 = coeff_c(:, 1);

        % 统一符号
        if dot(d0, e1) < 0
            d0 = -d0;
        end
        figure
        % ---------- 子图 1：三维视图 ----------
        subplot(1,2,1);
        scatter3(rel(:,1), rel(:,2), rel(:,3), 5, [0.7 0.7 0.7], 'filled');
        hold on;
        plot3(0, 0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

        % ★ PCA 三个主方向
        quiver3(0,0,0, e1(1)*L, e1(2)*L, e1(3)*L, ...
                'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        quiver3(0,0,0, e2(1)*L, e2(2)*L, e2(3)*L, ...
                'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        quiver3(0,0,0, n_local(1)*L, n_local(2)*L, n_local(3)*L, ...
                'm', 'LineWidth', 2, 'MaxHeadSize', 0.5);

        % 中心格子方向 d0（红）
        quiver3(0,0,0, d0(1)*L, d0(2)*L, d0(3)*L, ...
                'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);

        % 周围格子方向
        offsets = [-L, 0, L];
        for a = 1:3
            for b = 1:3
                if a == 2 && b == 2, continue; end
                c1 = offsets(a);
                c2 = offsets(b);
                in_cell = find(proj1 > c1 - L/2 & proj1 < c1 + L/2 & ...
                               proj2 > c2 - L/2 & proj2 < c2 + L/2);
                if length(in_cell) < 3, continue; end
                [coeff_j, ~, ~] = pca(rel(in_cell, :));
                d_j = coeff_j(:, 1);
                if dot(d_j, d0) < 0, d_j = -d_j; end

                cell_center_3d = c1 * e1 + c2 * e2;

                cos_theta = abs(dot(d_j, d0));
                cos_theta = min(max(cos_theta, -1), 1);
                theta_j = acosd(cos_theta);
                if theta_j < thetaThr
                    color = [0 0.7 0];
                else
                    color = [0.5 0.5 0.5];
                end

                quiver3(cell_center_3d(1), cell_center_3d(2), cell_center_3d(3), ...
                        d_j(1)*L/2, d_j(2)*L/2, d_j(3)*L/2, ...
                        'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
            end
        end
        axis equal;
        xlabel('X'); ylabel('Y'); zlabel('Z');
        title(sprintf('三维视图：点 %d, S_i=%d, n_{coh}=%d', ...
                      i, S(i), round(S(i)*8)));
        legend({'邻域点','中心点','e_1','e_2','n','d_0'}, 'Location', 'best');
        grid on;
        view(45, 30);

        % ---------- 子图 2：切平面视图 ----------
        subplot(1,2,2);
        scatter(proj1, proj2, 5, [0.7 0.7 0.7], 'filled');
        hold on;
        plot(0, 0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

        % ★ e1, e2 方向（就是坐标轴方向）
        quiver(0, 0, L, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        quiver(0, 0, 0, L, 'g', 'LineWidth', 2, 'MaxHeadSize', 0.5);

        % 九宫格边界
        edges = [-1.5*L, -0.5*L, 0.5*L, 1.5*L];
        for a = 1:4
            plot([edges(a) edges(a)], [-1.5*L 1.5*L], 'k-');
            plot([-1.5*L 1.5*L], [edges(a) edges(a)], 'k-');
        end

        % d0 投影（红）
        d0_proj = d0 - (d0' * n_local) * n_local;
        if norm(d0_proj) > eps
            d0_proj = d0_proj / norm(d0_proj);
            d0_proj_1 = d0_proj' * e1;
            d0_proj_2 = d0_proj' * e2;
            quiver(0, 0, d0_proj_1*L, d0_proj_2*L, ...
                   'r', 'LineWidth', 2, 'MaxHeadSize', 0.5);
        end

        % 周围格子方向（投影后）
        for a = 1:3
            for b = 1:3
                if a == 2 && b == 2, continue; end
                c1 = offsets(a);
                c2 = offsets(b);
                in_cell = find(proj1 > c1 - L/2 & proj1 < c1 + L/2 & ...
                               proj2 > c2 - L/2 & proj2 < c2 + L/2);
                if length(in_cell) < 3, continue; end
                [coeff_j, ~, ~] = pca(rel(in_cell, :));
                d_j = coeff_j(:, 1);
                if dot(d_j, d0) < 0, d_j = -d_j; end

                d_j_proj = d_j - (d_j' * n_local) * n_local;
                if norm(d_j_proj) < eps, continue; end
                d_j_proj = d_j_proj / norm(d_j_proj);
                d_j_1 = d_j_proj' * e1;
                d_j_2 = d_j_proj' * e2;

                cos_theta = abs(dot(d_j, d0));
                cos_theta = min(max(cos_theta, -1), 1);
                theta_j = acosd(cos_theta);
                if theta_j < thetaThr
                    color = [0 0.7 0];
                else
                    color = [0.5 0.5 0.5];
                end

                quiver(c1, c2, d_j_1*L/2, d_j_2*L/2, ...
                       'Color', color, 'LineWidth', 1.5, 'MaxHeadSize', 0.5);
            end
        end
        axis equal;
        xlabel('proj_1'); ylabel('proj_2');
        title('切平面视图：绿=一致，灰=不一致');
        legend({'邻域点','中心点','e_1','e_2','d_0'}, 'Location', 'best');
        grid on;
    end
end