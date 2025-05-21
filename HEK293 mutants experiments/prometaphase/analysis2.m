

%%
ff = figure(301);
ff.Position = [100, 100, 350, 350];  % 전체 그림 크기 조정

% 첫 번째 subplot
PCC_matrix = [PCC_wt(:,1), PCC_double(:,1), PCC_s262d(:,1)];
groupVector = repmat(1:size(PCC_matrix, 2), size(PCC_matrix, 1), 1);
groupVector = groupVector(:);
dataMatrix = PCC_matrix(:);
groupNames = {'Wild Type','T231D/S235D','S262D'};
plotSpread(PCC_matrix, 'distributionColors', 'k', 'markerSize', 7);
hold on;
boxplot(dataMatrix, groupVector, 'positions', [1 2 3], 'Widths', 0.6, 'Colors', 'k', 'Symbol', 'k+', 'Whisker', 1);
lines = findobj(gca, 'type', 'line');
set(lines, 'LineWidth', 1.5);
set(findobj(gca, 'tag', 'Median'), 'Color', 'red', 'LineWidth', 2);
ylabel('PCC value');
title('Tubulin-Tau')
set(gca, 'XTickLabel', groupNames);
ylim([0 1.1]);
xlim([0 4]);
hold off;

ff = figure(302);
ff.Position = [100, 100, 350, 350];  % 전체 그림 크기 조정

% 두 번째 subplot
PCC_matrix = [PCC_wt(:,2), PCC_double(:,2), PCC_s262d(:,2)];
groupVector = repmat(1:size(PCC_matrix, 2), size(PCC_matrix, 1), 1);
groupVector = groupVector(:);
dataMatrix = PCC_matrix(:);
plotSpread(PCC_matrix, 'distributionColors', 'k', 'markerSize', 7);
hold on;
boxplot(dataMatrix, groupVector, 'positions', [1 2 3], 'Widths', 0.6, 'Colors', 'k', 'Symbol', 'k+', 'Whisker', 1);
lines = findobj(gca, 'type', 'line');
set(lines, 'LineWidth', 1.5);
set(findobj(gca, 'tag', 'Median'), 'Color', 'red', 'LineWidth', 2);
ylabel('PCC value');
title('Centromere-Tau')
set(gca, 'XTickLabel', groupNames);
ylim([0 1.1]);
xlim([0 4]);
hold off;

