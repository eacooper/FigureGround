function resultsTable = create_eccentricity_comparison_table(stats, comparisons, meanDiffLabel, data, groups)
    % Create a formatted eccentricity comparison table with Cohen's D.
    % Inputs:
    %   stats        - the output from anova1
    %   comparisons  - the output from multcompare
    %   meanDiffLabel - the label for the mean difference column
    %   data         - the original data matrix (rows = observations, columns = groups)
    %   groups       - the group labels corresponding to the data

    % Extract original group names from stats without formatting
    originalGroupNames = stats.gnames;

    % Extract group indices for Ecc1 and Ecc2 from the comparisons
    Ecc1 = originalGroupNames(comparisons(:, 1));   % Group 1 Names
    Ecc2 = originalGroupNames(comparisons(:, 2));   % Group 2 Names

    % Format mean differences and CIs with 2 decimal places
    meanDiff = arrayfun(@(x) sprintf('%.2f', x), comparisons(:, 4), 'UniformOutput', false);
    ciLower = arrayfun(@(x) sprintf('%.2f', x), comparisons(:, 3), 'UniformOutput', false);
    ciUpper = arrayfun(@(x) sprintf('%.2f', x), comparisons(:, 5), 'UniformOutput', false);

    % Combine lower and upper CIs into a single string column
    CI = strcat("(", string(ciLower), ", ", string(ciUpper), ")");

    % Format p-values: Replace p < 0.001 with '<0.001'
    pValues = comparisons(:, 6);
    formattedPValues = strings(size(pValues));
    for i = 1:length(pValues)
        if pValues(i) < 0.001
            formattedPValues(i) = '<0.001';
        else
            formattedPValues(i) = sprintf('%.4f', pValues(i));
        end
    end

    % Compute Cohen's D for each comparison
    cohenD = zeros(size(comparisons, 1), 1);
    for i = 1:size(comparisons, 1)
        % Get group indices from the comparison
        group1_idx = comparisons(i, 1);
        group2_idx = comparisons(i, 2);

        % Extract data for the two groups being compared
        group1_data = data(strcmp(groups, originalGroupNames{group1_idx}));
        group2_data = data(strcmp(groups, originalGroupNames{group2_idx}));

        % Compute means and standard deviations
        M1 = mean(group1_data);
        M2 = mean(group2_data);
        SD1 = std(group1_data);
        SD2 = std(group2_data);
        n1 = length(group1_data);
        n2 = length(group2_data);

        % Pooled standard deviation
        pooledSD = sqrt(((n1 - 1) * SD1^2 + (n2 - 1) * SD2^2) / (n1 + n2 - 2));

        % Calculate Cohen's D and take the absolute value
        cohenD(i) = abs((M1 - M2) / pooledSD);
    end

    % Format Cohen's D to 2 decimal places
    cohenD = arrayfun(@(x) sprintf('%.2f', x), cohenD, 'UniformOutput', false);

    % Format the group names for display in the table
    formattedEcc1 = strrep(Ecc1, 'pt', '.');  % Replace "pt" with "."
    formattedEcc1 = strrep(formattedEcc1, '_deg', '');  % Remove "_deg"
    formattedEcc2 = strrep(Ecc2, 'pt', '.');  % Replace "pt" with "."
    formattedEcc2 = strrep(formattedEcc2, '_deg', '');  % Remove "_deg"

    % Create the results table with the appropriate column names
    resultsTable = table(formattedEcc1, formattedEcc2, meanDiff, CI, formattedPValues, cohenD, ...
        'VariableNames', {'Ecc1', 'Ecc2', meanDiffLabel, 'CI', 'p', 'CohensD'});
end
