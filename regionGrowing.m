function updatedPoints = regionGrowing(downsampledPtCloud, initialValues, distanceThreshold)
    % ptCloud - 点云数据 (N×3 矩阵)
    % initialValues - 每个点的初始赋值 (N×1 向量)
    % distanceThreshold - 距离阈值 (标定为0.1)

    % 确保 initialValues 是列向量
%     if isrow(initialValues)
%         initialValues = initialValues';
%     end
    
    % 获取点的总数
    numPoints = size(downsampledPtCloud.Location, 1);
    
    % 初始化点的标签
    labels = initialValues; % 标签：1、0.5、0
    
    % 查找所有赋值为1的点的索引
    seedIndices = find(labels == 1);
 
    % 创建一个空数组用于存储被更新的点的索引
    updatedIndices = [];
    
    % 匹配赋值为1的点与赋值为0.5的点
    while true
        changed = false;
        
        % 对于每个赋值为1的点
        for idx = seedIndices'
            % 创建一个半径为distanceThreshold的球形ROI
            roi = [downsampledPtCloud.Location(idx, 1) - distanceThreshold, downsampledPtCloud.Location(idx, 1) + distanceThreshold...
                downsampledPtCloud.Location(idx, 2) - distanceThreshold, downsampledPtCloud.Location(idx, 2) + distanceThreshold...
                downsampledPtCloud.Location(idx, 3) - distanceThreshold, downsampledPtCloud.Location(idx, 3) + distanceThreshold];
            
            % 查找在ROI内的点
            inROIIdx = findPointsInROI(downsampledPtCloud, roi);
            
            % 对于这些点，检查它们的赋值
            for j = inROIIdx'
                if labels(j) == 0.5
                    labels(j) = 1;
                    updatedIndices = [updatedIndices; j];
                    changed = true;
                end
            end
        end
        
        % 如果没有更多点需要修改，结束循环
        if ~changed
            break;
        end
        
        % 更新赋值为1的点的索引
        seedIndices = find(labels == 1);
    end
    
    % 返回更新后的点云数据
    updatedPoints = [downsampledPtCloud.Location, labels];
end
