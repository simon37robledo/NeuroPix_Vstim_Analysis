function [p_omnibus, Fobs, z_score] = omnibus_perm_test(respRates, baseRates, categories, nPerm)
if nargin < 4, nPerm = 5000; end

% Build groups: group 1 = baseline, groups 2..11 = categories 1..10
T = numel(respRates);
K = numel(unique(categories));  % should be 10

% Create vectors of values and group labels
% baseline is baseRates (T×1), stim values are respRates
vals = [baseRates(:); respRates(:)];          % length 2T
groups = [zeros(T,1); categories(:)];        % 0 = baseline, 1..K = categories
uniqueGroups = unique(groups);               % should be [0 1 .. K]

% Compute observed F-statistic (one-way ANOVA formula)
groupIDs = uniqueGroups;
gmean = arrayfun(@(g) mean(vals(groups==g)), groupIDs);
ng = arrayfun(@(g) sum(groups==g), groupIDs);
grandMean = mean(vals);
SSB = sum( ng .* (gmean - grandMean).^2 );
dfBetween = numel(groupIDs)-1;
SSW = sum( arrayfun(@(g) sum((vals(groups==g)-gmean(groupIDs==g)).^2), groupIDs) );
dfWithin = length(vals) - numel(groupIDs);
MSB = SSB/dfBetween;
MSW = SSW/dfWithin;
Fobs = MSB / MSW;

% Permutations
permF = zeros(nPerm,1);
n = length(vals);
for p = 1:nPerm
    idx = randperm(n);
    vals_perm = vals(idx);
    % reassign groups by taking same group sizes in original order
    % create a permuted grouping by shuffling labels
    groups_perm = groups(idx);
    % compute group means and sizes after permutation
    gmean_p = arrayfun(@(g) mean(vals_perm(groups_perm==g)), groupIDs);
    ng_p = arrayfun(@(g) sum(groups_perm==g), groupIDs);
    SSBp = sum( ng_p .* (gmean_p - mean(vals_perm)).^2 );
    SSWp = sum( arrayfun(@(g) sum((vals_perm(groups_perm==g)-gmean_p(groupIDs==g)).^2), groupIDs) );
    MSBp = SSBp/dfBetween;
    MSWp = SSWp/dfWithin;
    permF(p) = MSBp / MSWp;
end

p_omnibus = mean(permF >= Fobs);   % one-sided: larger F means groups differ
z_score = (Fobs - mean(permF)) / std(permF);

end
