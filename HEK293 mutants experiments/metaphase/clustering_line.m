% Step 1: Load Images
% Replace these paths with your actual file path
tau_file = 'Exp_11643_tau.tif';       % Tau channel
tubulin_file = 'Exp_11643_mt.tif';    % Tubulin (Spindle) channel
dna_file = 'Exp_11643_dna.tif';       % DNA channel

img_tau = imread(tau_file);
img_tubulin = imread(tubulin_file);
img_dna = imread(dna_file);

% Step 2: Define ROI and Crop Images
figure('Name', 'Select ROI on Tubulin Channel');
imshow(img_tubulin, []);
title('Draw a rectangle to define ROI. Double-click to finalize.');

% Draw ROI
roi = drawrectangle('Color', 'r', 'LineWidth', 2);
wait(roi);

% Get ROI position and crop region
roi_position = round(roi.Position); % [x, y, width, height]
x_start = max(1, roi_position(1));
y_start = max(1, roi_position(2));
x_end = min(size(img_tubulin, 2), x_start + roi_position(3) - 1);
y_end = min(size(img_tubulin, 1), y_start + roi_position(4) - 1);

cropped_tubulin = img_tubulin(y_start:y_end, x_start:x_end);
cropped_tau = img_tau(y_start:y_end, x_start:x_end);
cropped_dna = img_dna(y_start:y_end, x_start:x_end);

% Step 3: Region1 and Region2 Analysis
% Region1: Entire DNA region
dna_threshold = 0.2 * max(cropped_dna(:)); % Adjust threshold as needed
mask_region1 = cropped_dna > dna_threshold; % Mask for DNA regions
mask_region1 = bwconvhull(mask_region1, 'objects'); % Simplify mask (entire DNA region)

% Region2: Spindle region only
mt_threshold = 0.15 * max(cropped_tubulin(:)); % Threshold for spindle
mask_region2 = cropped_tubulin > mt_threshold; % Binary mask for spindle region
mask_region2_simplified = bwconvhull(mask_region2, 'objects'); % Simplify spindle region mask

scale_factor_y = 1.5; % Vertical scaling factor
[height, width] = size(mask_region2_simplified);

new_height = round(height * scale_factor_y);

mask_region2_resized_y = imresize(mask_region2_simplified, [new_height, width], 'nearest');

padded_mask_y = false(height, width); 

original_center_y = round(height / 2);
new_center_y = round(new_height / 2);

start_idx_y = max(1, original_center_y - new_center_y + 1);
end_idx_y = min(height, start_idx_y + new_height - 1);

if new_height > height
    cropped_resized_mask_y = mask_region2_resized_y(new_center_y - original_center_y + 1 : new_center_y + (height - original_center_y), :);
    padded_mask_y = cropped_resized_mask_y;
else
    padded_mask_y(start_idx_y:end_idx_y, :) = mask_region2_resized_y;
end
 
scale_factor_x = 1.1; % Horizontal scaling factor
[height_y, width_y] = size(padded_mask_y);

new_width = round(width_y * scale_factor_x);

mask_region2_resized_x = imresize(padded_mask_y, [height_y, new_width], 'nearest');

padded_mask_x = false(height_y, width_y); 

original_center_x = round(width_y / 2);
new_center_x = round(new_width / 2);

start_idx_x = max(1, original_center_x - new_center_x + 1);
end_idx_x = min(width_y, start_idx_x + new_width - 1);

if new_width > width_y
    cropped_resized_mask_x = mask_region2_resized_x(:, new_center_x - original_center_x + 1 : new_center_x + (width_y - original_center_x));
    padded_mask_x = cropped_resized_mask_x;
else
    padded_mask_x(:, start_idx_x:end_idx_x) = mask_region2_resized_x;
end

% final mask
mask_region2_resized_final = padded_mask_x;

% Step 5: Region1 excluding Region2
mask_region1_minus_region2 = mask_region1 & ~mask_region2_resized_final;

% DNA intensity values in each region
dna_region1 = cropped_dna(mask_region1);
dna_region2 = cropped_dna(mask_region2_resized_final);
dna_outside_region2 = cropped_dna(mask_region1_minus_region2);

dna_region1_minus_region2 = double(cropped_dna(mask_region1_minus_region2));
tau_region1_minus_region2 = double(cropped_tau(mask_region1_minus_region2));

if ~isempty(dna_region1_minus_region2) && ~isempty(tau_region1_minus_region2)
    % Normalize DNA and Tau intensities to the range [0, 1]
    dna_min = min(dna_region1_minus_region2);
    dna_max = max(dna_region1_minus_region2);
    tau_min = min(tau_region1_minus_region2);
    tau_max = max(tau_region1_minus_region2);

    dna_normalized = (dna_region1_minus_region2 - dna_min) / (dna_max - dna_min);
    tau_normalized = (tau_region1_minus_region2 - tau_min) / (tau_max - tau_min);

    % Scatter plot for normalized data
    fig_scatter = figure('Name', 'Scatter Plot: Draw Two Lines');
    scatter(dna_normalized, tau_normalized, 15, 'b', 'filled');
    hold on;
    xlabel('Normalized DNA Intensity');
    ylabel('Normalized Tau Intensity');
    title('Scatter Plot: Draw Two Lines');
    grid on;

    % Allow user to draw the first line
    disp('Draw the first line. Double-click to finish.');
    line1 = drawline('Color', 'r', 'LineWidth', 2);
    line1_pos = line1.Position; % [x1, y1; x2, y2]
    slope1 = (line1_pos(2, 2) - line1_pos(1, 2)) / (line1_pos(2, 1) - line1_pos(1, 1)); % Slope of line 1
    intercept1 = line1_pos(1, 2) - slope1 * line1_pos(1, 1); % Intercept of line 1

    % Allow user to draw the second line
    disp('Draw the second line. Double-click to finish.');
    line2 = drawline('Color', 'g', 'LineWidth', 2);
    line2_pos = line2.Position; % [x1, y1; x2, y2]
    slope2 = (line2_pos(2, 2) - line2_pos(1, 2)) / (line2_pos(2, 1) - line2_pos(1, 1)); % Slope of line 2
    intercept2 = line2_pos(1, 2) - slope2 * line2_pos(1, 1); % Intercept of line 2
    
    % Calculate intersection of the two lines
    A = [-slope1, 1; -slope2, 1];
    b = [intercept1; intercept2];
    intersection = A \ b; % [x_intersect; y_intersect]
    new_origin = intersection'; % Use the intersection point as the new origin
    disp(['New origin (intersection) calculated at: (', num2str(new_origin(1)), ', ', num2str(new_origin(2)), ')']);

    
    % Step 5: Generate Final Figure
    figure('Name', 'MT Threshold, Region Analysis, and Cross-Correlation', 'NumberTitle', 'off');

    % Original DNA Image
    subplot(2, 5, 1);
    imshow(cropped_dna, []);
    title('Original DNA Image');

    % Region1 Mask
    subplot(2, 5, 2);
    imshow(mask_region1, []);
    title('Region1: Entire DNA Region');

    % Original Microtubule Image
    subplot(2, 5, 3);
    imshow(cropped_tubulin, []);
    title('Original Microtubule Image');

    % Region2 Mask
    subplot(2, 5, 4);
    imshow(mask_region2_resized_final, []);
    title('Region2: Spindle Region');

    % Region1 Minus Region2 Mask
    subplot(2, 5, 5);
    imshow(mask_region1_minus_region2, []);
    title('Region1 Minus Region2');

    % Scatter Plot
    subplot(2, 5, 6);
    scatter(dna_normalized, tau_normalized, 10, 'b', 'filled');
    hold on;
    
    % Plot line1
    line1_pos = line1.Position; % Get the position of the first line
    plot(line1_pos(:, 1), line1_pos(:, 2), 'r-', 'LineWidth', 2); % Plot the first line in red
    
    % Plot line2
    line2_pos = line2.Position; % Get the position of the second line
    plot(line2_pos(:, 1), line2_pos(:, 2), 'g-', 'LineWidth', 2); % Plot the second line in green

    % Highlight selected point
    scatter(new_origin(1), new_origin(2), 50, 'k', 'filled'); 

    hold off;

    xlim([0,1]);
    ylim([0,1]);
    xlabel('Normalized DNA Intensity');
    ylabel('Normalized Tau Intensity');
    
    title('Scatter Plot with Selected New Origin');
    grid on;

    % Step 6: K-means Clustering Based on New Origin
    subplot(2, 5, 7);
    
    % Shift data relative to the new origin for cosine distance calculation only
    dna_shifted = dna_normalized - new_origin(1);
    tau_shifted = tau_normalized - new_origin(2);
    
    % Combine shifted data for cosine distance calculation
    shifted_data = [dna_shifted, tau_shifted];
    
    % Normalize shifted data to unit vectors for cosine similarity
    shifted_data_normalized = normalize(shifted_data, 2);
    
    % K-means clustering with cosine distance
    num_clusters = 2; % Number of clusters
    [cluster_indices, cluster_centers] = kmeans(shifted_data_normalized, num_clusters, ...
        'Distance', 'cosine', 'Replicates', 10);
    
    % Scatter plot for clustered data
    hold on;
    colors = lines(num_clusters);
    slopes = zeros(num_clusters, 1); % Store slopes for each cluster

    for cluster_id = 1:num_clusters
        cluster_points = (cluster_indices == cluster_id);
        cluster_dna = dna_normalized(cluster_points);
        cluster_tau = tau_normalized(cluster_points);
        
        % Scatter plot for the cluster
        scatter(cluster_dna, cluster_tau, 15, colors(cluster_id, :), 'filled');
        
        % Calculate trend line for the cluster
        if numel(cluster_dna) > 1
            p = polyfit(cluster_dna, cluster_tau, 1); % Linear fit
            slopes(cluster_id) = p(1); % Store slope
            x_fit = linspace(min(cluster_dna), max(cluster_dna), 100);
            y_fit = polyval(p, x_fit);
            plot(x_fit, y_fit, '-', 'Color', 'k', 'LineWidth', 1.5);
        end
    end
    
    xlim([0,1]);
    ylim([0,1]);
    xlabel('Normalized DNA Intensity');
    ylabel('Normalized Tau Intensity');
    title('K-means Clustering with New Origin', 'FontSize', 10);
    grid on;
    hold off;

    % Identify the cluster with the larger slope
    [~, max_slope_cluster_id] = max(slopes); % Cluster with the maximum slope

    % Step 7: Cross-Correlation Map
    subplot(2, 5, 8);
    cc_map_region1_minus_region2 = zeros(size(cropped_tau));
    
    % Filter data points belonging to the selected cluster
    selected_cluster_points = (cluster_indices == max_slope_cluster_id);
    filtered_dna_cc = dna_region1_minus_region2(selected_cluster_points);
    filtered_tau_cc = tau_region1_minus_region2(selected_cluster_points);
    
    % Update CC map with filtered cluster data
    filtered_indices = false(size(mask_region1_minus_region2));
    filtered_indices(mask_region1_minus_region2) = selected_cluster_points;
    cc_map_region1_minus_region2(filtered_indices) = ...
        filtered_dna_cc .* filtered_tau_cc;
    
    imagesc(cc_map_region1_minus_region2, [0, max(cc_map_region1_minus_region2(:))]);
    colormap('hot');
    colorbar;
    title('DNA and Tau Cross-Correlation Map (Filtered)', 'FontSize', 10);
    xlabel('X (Pixels)');
    ylabel('Y (Pixels)');
    
    % Step 8: Heat Map
    subplot(2, 5, 9);
    
    % Define histogram range for the filtered cluster data
    xedges_cc = linspace(double(min(filtered_dna_cc)), double(max(filtered_dna_cc)), 50);
    yedges_cc = linspace(double(min(filtered_tau_cc)), double(max(filtered_tau_cc)), 50);
    
    % Compute 2D histogram for the filtered cluster data
    [counts_cc, xedges_cc, yedges_cc] = histcounts2(filtered_dna_cc, filtered_tau_cc, xedges_cc, yedges_cc);
    
    % Normalize counts for visualization
    counts_cc_sum_normalized = counts_cc / sum(counts_cc(:));
    
    % Compute heatmap centers
    xCenters_cc = (xedges_cc(1:end-1) + xedges_cc(2:end)) / 2;
    yCenters_cc = (yedges_cc(1:end-1) + yedges_cc(2:end)) / 2;
    
    % Plot heatmap
    imagesc(xCenters_cc, yCenters_cc, counts_cc_sum_normalized');
    axis xy;
    xlabel('DNA Intensity');
    ylabel('Tau Intensity');
    colorbar;
    title('Heat Map: DNA vs Tau Intensity');
    
else
    disp('No valid data points available.');
end



% Step 9: Proportion Calculations
% Calculate total sums for Region1 (Entire DNA Region)
total_dna_sum = sum(cropped_dna(mask_region1)); % Sum of DNA intensities in Region1
total_tau_sum = sum(cropped_tau(mask_region1)); % Sum of Tau intensities in Region1

% Update proportion calculations for filtered cluster
filtered_dna_sum = sum(filtered_dna_cc); % DNA intensity sum for selected cluster
filtered_tau_sum = sum(filtered_tau_cc); % Tau intensity sum for selected cluster


filtered_pixels_dna = numel(filtered_dna_cc); % Number of DNA data points in selected cluster
filtered_pixels_tau = numel(filtered_tau_cc); % Number of Tau data points in selected cluster

percentage_filtered_dna = (filtered_dna_sum / total_dna_sum) * 100;
percentage_filtered_tau = (filtered_tau_sum / total_tau_sum) * 100;

% Command Window 
disp('=== Filtered Cluster Proportion Results ===');
disp(['DNA Intensity Proportion: ', sprintf('%.2f', percentage_filtered_dna), '%']);
disp(['Tau Intensity Proportion: ', sprintf('%.2f', percentage_filtered_tau), '%']);

% Display Results on the Final Plot
annotation('textbox', [0.75, 0.1, 0.2, 0.3], 'String', {
    sprintf('DNA Intensity Proportion: %.2f%%', percentage_filtered_dna), ...
    sprintf('Tau Intensity Proportion: %.2f%%', percentage_filtered_tau)}, ...
    'FitBoxToText', 'on', 'BackgroundColor', 'white', 'FontSize', 10);
