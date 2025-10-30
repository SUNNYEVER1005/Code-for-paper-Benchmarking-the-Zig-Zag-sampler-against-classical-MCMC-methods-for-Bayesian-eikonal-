function burnInPeriod = burned_estimate(data, alpha)
% burned_estimate  Estimate burn-in period for each row of data and return the maximum
%
%   burnInPeriod = burned_estimate(data, alpha)
%
%   Inputs:
%     data  – 3×n numeric matrix, each row is a time series
%     alpha – tolerance threshold for mean change
%
%   Output:
%     burnInPeriod – the maximum burn-in index across the three rows

    [numRows, numCols] = size(data);
    if numRows ~= 3
        error('Input data must be a 3×n matrix.');
    end

    burnIns = zeros(1, numRows);

    % For each row, estimate its burn-in
    for iRow = 1:numRows
        rowData      = data(iRow, :);
        currentMean  = inf;
        proposalMean = 0;
        t            = 0;

        for t = 1:numCols
            % Every 100 points compute the windowed mean
            if mod(t, 100) == 0
                % Compute mean of the last 100 points
                windowStart    = t - 99;
                proposalMean   = mean(rowData(windowStart:t));
                tolerance      = abs(proposalMean - currentMean);

                % Stop if change is within alpha
                if tolerance <= alpha
                    break;
                end

                currentMean = proposalMean;
            end
        end

        burnIns(iRow) = t;
    end

    % Return the maximum burn-in across the three rows
    burnInPeriod = max(burnIns);
end
