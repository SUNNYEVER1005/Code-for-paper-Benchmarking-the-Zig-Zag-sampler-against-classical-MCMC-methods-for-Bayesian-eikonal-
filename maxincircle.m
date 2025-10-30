function Tmax = maxincircle(interpFunc, pos, r)
    sampleDensity = 10;
    x0 = pos(1); z0 = pos(2);
    xs = linspace(x0 - r, x0 + r, sampleDensity);
    zs = linspace(z0 - r, z0 + r, sampleDensity);
    [X, Z] = ndgrid(xs, zs);
    mask = (X - x0).^2 + (Z - z0).^2 <= r^2;
    
    Tmax = cellfun(@(f) getTmax(f, X, Z, mask), interpFunc);
end

function tmax = getTmax(f, X, Z, mask)
    % Compute the maximum absolute value of the given interpolant \( f \) 
    % inside the circle defined by the X–Z grid.
    tvals = f(X, Z);
    tmax = max(abs(tvals(mask)));
end  








% function Tmax = maxincircle(i,N_r,interpFunc, pos, r)
% % interpMaxInCircle 在以 (x0,z0) 为中心半径为 r 的圆内查找插值函数 interpFunc 的最大值
% %
% % 输入参数：
% %   i             : 输入插值函数或矩阵的第i个指标
% %   interpFunc    : 已构造的 griddedInterpolant 插值函数
% %   pos=[x0, z0,\tau]        : 查询点坐标
% %   r             : 圆的半径
% %   sampleDensity : (可选) 圆内采样点的密度（默认值：50）
% %
% % 输出参数：
% %   Tmax          : 圆内的最大 T 值
% 
%     % if nargin < 5
%     %     sampleDensity = 50; % 默认采样密度
%     % end
%     sampleDensity = 10; % 默认采样密度
% 
%     % 在查询点周围构造一个包围正方形的网格
%     x0=pos(1);z0=pos(2);  %查询点
% 
%     xs = linspace(x0 - r, x0 + r, sampleDensity);
%     zs = linspace(z0 - r, z0 + r, sampleDensity);
%     [X, Z] = ndgrid(xs, zs);
% 
%     % 找出落在圆内的采样点
%     mask = (X - x0).^2 + (Z - z0).^2 <= r^2;
% 
%     % 对所有网格点进行插值计算
%     Tvals = interpFunc{i}(X, Z);
% 
%     % 提取圆内的值并计算最大值
%     Tmax = max(abs(Tvals(mask)));
% end