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
