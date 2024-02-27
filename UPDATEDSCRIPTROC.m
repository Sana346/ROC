% Clear all variables, close all figures, and clear the command window
clear all;
close all;
clc;

% Display progress message
disp('PROGRESS: Script initialized');

%% IMPORT EXTERNAL DATA
image = imread('TI AND TiW_30KEV_300X_006 bs.tif');
disp('PROGRESS: Image imported');

%% CONVERT INPUT IMAGE TO GRAYSCALE, CROP IT, AND RESIZE IT
image_gray = im2gray(image);
image_crop = imcrop(image_gray, [0.5 0.5 1023.49803921569 888.629411764706]);
J = image_crop;
histogramValues = imhist(image_crop);

%% BINARIZE PROCESSED IMAGE AND EXTRACT & RENDER ITS BOUNDARIES
BW = imbinarize(J,'adaptive','ForegroundPolarity','dark','Sensitivity',0.03);
[B, ~] = bwboundaries(BW,'holes');

min_radii = cell(1, numel(B));
max_radii = cell(1, numel(B));
radii_for_hist = cell(1, numel(B));
num_ignored_pores = 0;

% Open the file for writing
fileID = fopen('pore_coordinates.txt','w');

% Check if the file was successfully opened
if fileID == -1
    error('Unable to open the file for writing.');
end

for k = 1:length(B)
    boundary = B{k};
    x_values = boundary(:, 2);
    y_values = boundary(:, 1);
end
    % Determine if this particular pore boundary should be ignored
    if any(x_values <= 1) || any(x_values >= size(image_crop, 2)) || any(y_values <= 1) || any(y_values >= size(image_crop, 1))
        % Skip this boundary as it's on the border
        num_ignored_pores = num_ignored_pores + 1;
%         continue;

    end
%       plot(x_values, y_values, 'r', 'LineWidth', 2);
%     hold off;
%     title(sprintf('Pore boundary %d', k));
%     xlabel(['\it{x}', '\rm', sprintf('-coordinate of Pore %d (px)', k)], 'FontSize', 12);
%     ylabel(['\it{y}', '\rm', sprintf('-coordinate of Pore %d (px)', k)], 'FontSize', 12)
    % Write pore coordinates to the text file
    fprintf(fileID,'Pore %d:\n', k);
    fprintf(fileID,'x values: %s\n', num2str(x_values'));
    fprintf(fileID,'y values: %s\n', num2str(y_values'));
    fprintf(fileID,'\n');

    %% Calculate curvature manually
    dx = diff(x_values);
    dy = diff(y_values);
    ds = sqrt(dx.^2 + dy.^2);
    curvature = 2 * abs(dx(2:end) .* dy(1:end-1) - dx(1:end-1) .* dy(2:end)) ./ (ds(1:end-1) .* ds(2:end));

    % Calculate the radii of curvature & key statistics
    radii_for_hist{k} = 1 ./ curvature;

    % Find the minimum curvature value excluding the endpoints
    curvature_no_endpoints = curvature(2:end-1);
    if ~isempty(curvature_no_endpoints)
        min_radii{k} = 1 / min(curvature_no_endpoints);
    else
        min_radii{k} = NaN; % Handle the case when there are no valid curvature values
    end

    % Find the maximum curvature value, if the curvature vector is not empty
    if ~isempty(curvature)
        max_radii{k} = 1 / max(curvature);
    else
        max_radii{k} = NaN; % Handle the case when there are no valid curvature values
    end

    fprintf('Min & max radii of curvature for Pore %d are %f & %f\n', k, min_radii{k}, max_radii{k});
    
% %     %% Create a histogram of radii_for_hist{k}     

    figure;
%     subplot(1, 2, 1);
%     plot(radii_for_hist{k},'-o');
%     title(['Radii of curvature for pore ', num2str(k)]);
%     xlabel('Data point index');
%     ylabel('Radius of curvature');
%     grid on;
 subplot(1, 2, 2);
    histogram(radii_for_hist{k}, 20); 
    title(['Histogram of radii of curvature for Pore ', num2str(k)]);
    xlabel('Radius of curvature');
    ylabel('Number frequency');
    grid on;
    set(gcf, 'Position', get(0, 'Screensize'));
%     Overlay pore boundaries onto the original image
%    % figure;
%% Overlay pore boundaries onto the original image
% %     figure;
    subplot(1,2,1)
    imshow(image_crop);
    hold on;
    plot(x_values, y_values, 'r', 'LineWidth', 2);
    hold off;
    title(sprintf('Pore boundary %d', k));
    xlabel(['\it{x}', '\rm', sprintf('-coordinate of pore %d (px)', k)], 'FontSize', 12);
    ylabel(['\it{y}', '\rm', sprintf('-coordinate of pore %d (px)', k)], 'FontSize', 12);
    pause;
 
%     imshow(image_crop);
% %     hold on;
% % for k = 1:length(B)
% %     boundary = B{k};
% %     x_values = boundary(:, 2);
% %     y_values = boundary(:, 1);
% %     if any(x_values <= 1) || any(x_values >= size(image_crop, 2)) || any(y_values <= 1) || any(y_values >= size(image_crop, 1))
% %         continue;
% %     end
% %     hold on;
% %     plot(x_values, y_values, 'r', 'LineWidth', 2);
% % end
%     hold off;
%     title(sprintf('Pore boundary %d', k));
%     xlabel(['\it{x}', '\rm', sprintf('-coordinate of Pore %d (px)', k)], 'FontSize', 12);
%     ylabel(['\it{y}', '\rm', sprintf('-coordinate of Pore %d (px)', k)], 'FontSize', 12)  
% %     Create a subplot for the current pore's boundary on the original image
% %         subplot(length(B), 2, (k-1)*2+2);
% %         imshow(image_crop);
% %         hold on;
% %         plot(x_values, y_values, 'r', 'LineWidth', 2);
% %         hold off;
% %         title(['Pore boundary ', num2str(k)]);
% %         xlabel(['\it{x}', '\rm', sprintf('-coordinate of Pore %d (px)', k)], 'FontSize', 12);
% %         ylabel(['\it{y}', '\rm', sprintf('-coordinate of Pore %d (px)', k)], 'FontSize', 12);
% 
% %     end
% % end
% %     subplot(1, 2, 2);
% %     histogram(radii_for_hist{k}, 20); 
% %     title(['Histogram of radii of curvature for Pore ', num2str(k)]);
% %     xlabel('Radius of curvature');
% %     ylabel('Number frequency');
% %     grid on;
% %     set(gcf, 'Position', get(0, 'Screensize'));
% %     
% %     %% Overlay pore boundaries onto the original image
% %     figure;
% %     imshow(image_crop);
% %     hold on;
% %     plot(x_values, y_values, 'r', 'LineWidth', 2);
% %     hold off;
% %     title(sprintf('Pore boundary %d', k));
% %     xlabel(['\it{x}', '\rm', sprintf('-coordinate of pore %d (px)', k)], 'FontSize', 12);
% %     ylabel(['\it{y}', '\rm', sprintf('-coordinate of pore %d (px)', k)], 'FontSize', 12);
%     pause;
%     
% % end

fclose(fileID);

% Display the number of pores ignored as they intersect with the image boundary
disp(['Number of pores ignored as they intersect with the image boundary = ', num2str(num_ignored_pores), ' out of ', num2str(length(B))]);

% Calculate the length and width of the binary image
BW_length = size(BW, 1);
BW_width = size(BW, 2);

%% POROSITY
% Calculate porosity using the formula
porosity = (BW_length * BW_width - sum(BW(:))) / (BW_length * BW_width);
fprintf('Porosity: %f%%\n', porosity * 100);

%% 
%% MASTER FIGURE

% Create master figure & panels
figure;
subplot(2, 2, 1);
imshow(image_crop);
title('1. Original grayscale image');

subplot(2, 2, 2);
bar(histogramValues);
title('2. Brightness histogram');
xlabel('Pixel value');
ylabel('Pixel frequency');

subplot(2, 2, 3);
imshow(BW);
title('3. Binarised image');

subplot(2, 2, 4);
imshow(image_crop);
hold on;
for k = 1:length(B)
    boundary = B{k};
    x_values = boundary(:, 2);
    y_values = boundary(:, 1);
    if any(x_values <= 1) || any(x_values >= size(image_crop, 2)) || any(y_values <= 1) || any(y_values >= size(image_crop, 1))
        continue;
    end
    plot(x_values, y_values, 'r', 'LineWidth', 2);
end
hold off;
title('4. Pore boundaries on original image');
set(gcf, 'Position', get(0, 'Screensize'));

% Display completion message
disp('PROGRAM SUCCESSFULLY COMPLETED');
