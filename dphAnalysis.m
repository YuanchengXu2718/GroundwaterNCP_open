% Author: Yuancheng Xu
% Affiliation: Prof. Di Long's group at Tsinghua University, Beijing, PRC
% Email: xuyc24@mails.tsinghua.edu.cn
% Last updated: 25/04/2025

% Code for spatial interpolation and uncertainty analysis of regionally 
% average groundwater depth anomaly

tic; % record the time used for the script

% parameters
nSample = 0.95; % percentage of sample chosen when doing Monte Carlo experiments
maxLoop = 10000; % nLoop times of Monte Carlo experiments
ci = 0.9; % confidence interval
% 20 years of data from Jan. 2005 to Dec. 2024
timeseries = datetime(2005, 1, 1):calmonths(1):datetime(2024, 12, 1);
timeseries = timeseries';

% parallel computing setting    
CoreNum = 24; % set the number of cores used in parallel computing
if isempty(gcp('nocreate'))
    parpool(CoreNum);
end

f = waitbar(0, 'Loading...');
root = "E:\2024groundwater_insitu\In_situ_writing_revise1\openData";
% read the Excel file of groundwater depth at wells
path_input = root + "\GroundwaterDepth.xlsx";
raw = readmatrix(path_input);

% type of monitoring well
% if 10 < type <= 19, the well monitors unconfined aquifers
% type 11 represents the first unconfined aquifers
% type 12 represents the second unconfined aquifers, etc.
% if 20 < type <= 29, the well monitors confined aquifers
% type 21 represents the first confined aquifers
% type 22 represents the second confined aquifers, etc.
id = raw(:, 1);
type = raw(:, 3);
% longitude and latitude of the monitoring well
lon = raw(:, 4);
lat = raw(:, 5);
% monthly depth to groundwater depth data
dphdata = raw(:, 6:end);

% change in depth between two consecutive months
data = zeros(size(dphdata));

% preprocessing parameters
intrplGap = 3; % maximum allowable gap for interpolation
minLen = 48; % minimum amount of data required
windowOL = 24; % moving window length for identifying outliers


% type of well selected
% when doing interpolation, we only include wells in the first layer of
% unconfined (type = 11) and confined (type = 21) aquifers
typeMask = type == 21;

% load the mask for interpolation
% choose your region of interest here
path_map = root + "\mask\ncp.tif"; % North China Plain
% path_map = root + "\mask\city_beijing.tif"; % Beijing
% path_map = root + "\mask\city_baoding.tif"; % Baoding
% path_map = root + "\mask\city_hengshui.tif"; % Hengshui
% path_map = root + "\mask\city_tangshan.tif"; % Tangshan
% property of the mask
[map, R] = readgeoraster(path_map);

% -------------------------------------------------------------------------
% Start of preprocessing

% temporal interpolation (linear)
for i = 1:size(dphdata, 1)
    temInter = 0;
    for j = 1:size(dphdata, 2)
        if temInter > 0
            if ~isnan(dphdata(i, j)) && j == temInter + 1 % valid, valid
                temInter = j;
            elseif ~isnan(dphdata(i, j)) && j <= temInter + intrplGap + 1 % valid, nan (<= intrplGap), valid 
                dphdata(i, temInter:j) = linspace(dphdata(i, temInter), dphdata(i, j), j - temInter + 1);
                temInter = j;
            elseif ~isnan(dphdata(i, j)) && j > temInter + intrplGap + 1 % valid, nan (> intrplGap), valid
                temInter = j;
            end
        elseif ~isnan(dphdata(i, j))
            temInter = j;
        else
            continue
        end
    end
end

for i = 1:size(dphdata, 1)
    % mask out wells without enough data
    if sum(~isnan(dphdata(i, :))) < minLen
        dphdata(i, :) = nan;
    end
    % outlier detection
    outlierMask = zeros(size(dphdata, 2), 1);
    for j = windowOL:size(dphdata, 2)
        dph_window = dphdata(i, j - windowOL + 1 : j);
        avg_window = mean(dph_window, 'omitmissing');
        sigma_window = sqrt(var(dph_window, 'omitmissing'));
        mask_window = dph_window > avg_window + 3*sigma_window |...
                        dph_window < avg_window - 3*sigma_window;
        outlierMask(j - windowOL + 1 : j) = mask_window' + outlierMask(j - windowOL + 1 : j);
    end
    outlierMask = outlierMask > 0.5 * windowOL;
    dphdata(i, outlierMask) = nan;
end

% End of preprocessing
% -------------------------------------------------------------------------

% delta depth at each well
for i = 2:size(data, 2)
    data(:, i) = dphdata(:, i) - dphdata(:, i-1);
end

% reproject the locations of wells to UTM zone 50N (EPSG: 32650)
projUTM = projcrs(32650);
[utm_x, utm_y] = projfwd(projUTM, lat, lon);

% coordinates of the entire map
map_x = linspace(R.XWorldLimits(1), R.XWorldLimits(2), size(map, 2));
map_y = linspace(R.YWorldLimits(2), R.YWorldLimits(1), size(map, 1));
map_res = R.CellExtentInWorldX;

% find the coordinates of pixels in the selected region
nPixel = sum(map, 'all');
map_xoi = zeros(nPixel, 1);
map_yoi = zeros(nPixel, 1);
loc_map = zeros(nPixel, 2);
count = 0;
for i = 1:size(map, 1)
    for j = 1:size(map, 2)
        if map(i, j) == 1
            count = count + 1;
            map_xoi(count) = map_x(j);
            map_yoi(count) = map_y(i);
            loc_map(count, :) = [i, j];
        end
    end
end

% find the location of in-situ groundwater depth in the map
waitbar(0, f, 'Start interpolation...');
deltaDph = zeros(size(timeseries, 1), 1);
cumDph = zeros(size(timeseries, 1), 1);

% parfor is used for parallel computing
parfor step = 2:size(timeseries, 1)
    % find available wells in the given month
    wellMask = typeMask & ~isnan(data(:, step));
    well_oi = data(wellMask, step);
    wellX_oi = utm_x(wellMask);
    wellY_oi = utm_y(wellMask);
    
    % spatial interpolation using Thiessen polygons (nearest neighbours)
    mapNN = GW_Dph_NN(well_oi, wellX_oi, wellY_oi, map_x, map_y, map_res, map, loc_map);

    % calculate the weighted average depth
    deltaDph(step) = mean(mapNN(mapNN~=-99), 'all', 'omitmissing');
end

for step = 2:size(timeseries, 1)
    cumDph(step) = cumDph(step - 1) + deltaDph(step);
end

% Uncertainty analysis

% record results of Monte Carlo simulation
deltaDph_loop = zeros(size(timeseries, 1), maxLoop);
cumDph_loop = zeros(size(timeseries, 1), maxLoop);

% Except for randomly selecting nSample of the input wells,
% everything else is the same as the code before
for nLoop = 1:maxLoop
    waitbar(nLoop/maxLoop, f, ['nLoop = ', num2str(nLoop)]);

    % randomly select nSample of the input wells
    dphdataR = dphdata;
    for step = 1:size(timeseries, 1)
        randomMask = rand(size(dphdata, 1), 1) > nSample;
        dphdataR(randomMask, step) = nan;
    end
    % change in depth between two consecutive months
    dataR = zeros(size(dphdata));
    for i = 2:size(dataR, 2)
        dataR(:, i) = dphdataR(:, i) - dphdataR(:, i-1);
    end

    parfor step = 1:size(timeseries, 1)
        wellMask = typeMask & ~isnan(dataR(:, step));
        well_oi = data(wellMask, step);
        wellX_oi = utm_x(wellMask);
        wellY_oi = utm_y(wellMask);

        mapNN = GW_Dph_NN(well_oi, wellX_oi, wellY_oi, map_x, map_y, map_res, map, loc_map);

        deltaDph_loop(step, nLoop) = mean(mapNN(mapNN~=-99), 'all', 'omitmissing');
    end
end
close(f);
for step = 2:size(timeseries, 1)
    cumDph_loop(step, :) = cumDph_loop(step - 1, :) + deltaDph_loop(step, :);
end

% % calculate the upper and lower bound of confidence interval
sorted = sort(cumDph_loop, 2, 'ascend');
lowerLoop = sorted(:, floor(0.5*(1-ci)*maxLoop)+1 );
upperLoop = sorted(:, ceil(0.5*(1+ci)*maxLoop) );

% visualize the results
close all;
figure(1);
plot(timeseries, cumDph - mean(cumDph), 'r-', 'LineWidth', 1);
hold on;
fill([timeseries; flipud(timeseries)], [lowerLoop - mean(cumDph); flipud(upperLoop) - mean(cumDph)],...
    'r', 'FaceAlpha', 0.3, 'EdgeColor','none');
hold on;
legend('Groundwater depth anomaly', 'Uncertainty (Thiessen)');
ylabel('Groundwater depth anomaly (m)');
set (gca,'YDir','reverse');
hold off;

% find the final result here:
% regionally average depth, lower bound of confidence interval, and
% upper bound of confidence interval
AttentionForResult = [cumDph - mean(cumDph), lowerLoop - mean(cumDph), upperLoop - mean(cumDph)];
AttentionForResult_cum = [cumDph, lowerLoop, upperLoop];

toc; % record the time used for the script
