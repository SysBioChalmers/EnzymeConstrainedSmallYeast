function proteinContent = linearProteinContent(growthRate)
    %Load protein content data
    x = [0.025 0.05 0.10 0.15 0.20 0.25 0.28 0.3 0.35 0.40];
    y = 0.01 * [37.9228 40.0524 42.0712 42.827 46.8089 49.2455 53.3415 47.3331 46.9645 45.7556];

    x = [0.025 0.28 0.4];
    y = [0.379 0.533 0.457];
    
    proteinContent = interp1q(x', y', growthRate');
    proteinContent(growthRate < min(x)) = y(1);
    proteinContent(growthRate > max(x)) = y(end);    
    
end