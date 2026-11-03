function thr = otsuThreshold(x)
% Otsu 阈值法（调用时先排除 0 值）

    x = x(~isnan(x));
    if isempty(x)
        thr = 0;
        return;
    end

    x_min = min(x);
    x_max = max(x);
    if x_max == x_min
        thr = x_min;
        return;
    end
    x_norm = (x - x_min) / (x_max - x_min);

    nbins = 256;
    [counts, centers] = hist(x_norm, nbins);
    p = counts / sum(counts);

    mu_total = sum(centers .* p);
    omega = cumsum(p);
    mu = cumsum(centers .* p);
    sigma_b_sq = (mu_total * omega - mu).^2 ./ (omega .* (1 - omega) + eps);
    [~, idx] = max(sigma_b_sq);
    thr_norm = centers(idx);

    thr = thr_norm * (x_max - x_min) + x_min;
end