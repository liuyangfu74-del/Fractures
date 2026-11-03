%% 均值位移聚类函数（无过滤版本）
function [labels, cluster_centers] = meanShift(data, bandwidth)
    % 基于文献描述的均值位移聚类算法
    % 输入:
    %   data - 输入数据 (n x 3), 法向量数据
    %   bandwidth - 带宽参数
    % 输出:
    %   labels - 聚类标签
    %   cluster_centers - 聚类中心
    
    [n, dim] = size(data);
    labels = zeros(n, 1);
    cluster_centers = [];
    cluster_count = 0;
    
    fprintf('开始均值位移聚类，数据点数量: %d, 带宽: %.3f\n', n, bandwidth);
    
    % 归一化数据（法向量已经在单位球面上）
    data_norm = data ./ vecnorm(data, 2, 2);
    
    % 初始化访问标记
    visited = false(n, 1);
    convergence_threshold = 1e-5 * bandwidth;
    
    % 进度显示
    fprintf('聚类进度: ');
    
    for i = 1:n
        if ~visited(i)
            % 显示进度
            if mod(i, max(1, floor(n/10))) == 0
                fprintf('.');
            end
            
            % 以当前点作为初始中心
            current_center = data_norm(i, :);
            iteration = 0;
            max_iterations = 100;
            
            % 均值位移迭代
            while iteration < max_iterations
                iteration = iteration + 1;
                
                % 计算当前中心与所有点的角度距离（球面距离）
                cos_distances = data_norm * current_center';
                angular_distances = acos(min(max(cos_distances, -1), 1)); % 防止数值误差
                
                % 找到在带宽范围内的点（使用角度距离）
                in_bandwidth = angular_distances <= bandwidth;
                
                if sum(in_bandwidth) == 0
                    break;
                end
                
                % 计算新的中心（球面均值）
                points_in_window = data_norm(in_bandwidth, :);
                new_center = mean(points_in_window, 1);
                new_center = new_center / norm(new_center); % 归一化到单位球面
                
                % 检查收敛（使用角度变化）
                center_shift = acos(min(max(dot(current_center, new_center), -1), 1));
                
                if center_shift < convergence_threshold
                    break;
                end
                
                current_center = new_center;
            end
            
            % 检查是否与现有聚类中心合并
            merged = false;
            merge_threshold = bandwidth / 2;
            
            for k = 1:cluster_count
                % 计算与现有聚类中心的角度距离
                center_distance = acos(min(max(dot(current_center, cluster_centers(k, :)), -1), 1));
                
                if center_distance < merge_threshold
                    % 合并到现有聚类
                    labels(in_bandwidth) = k;
                    visited(in_bandwidth) = true;
                    merged = true;
                    break;
                end
            end
            
            if ~merged
                % 创建新聚类
                cluster_count = cluster_count + 1;
                cluster_centers(cluster_count, :) = current_center;
                labels(in_bandwidth) = cluster_count;
                visited(in_bandwidth) = true;
            end
        end
    end
    
    fprintf('\n聚类完成，共发现 %d 个聚类\n', cluster_count);
end