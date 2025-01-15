% Author: Yuancheng Xu
% Affiliation: Prof. Di Long's group at Tsinghua University, Beijing, PRC
% Email: xuyc24@mails.tsinghua.edu.cn
% Last updated: 15/01/2025

% Code for spatial interpolation and uncertainty analysis of regionally 
% average groundwater depth

tic; % record the time used for the script

% parameters
nSample = 0.9; % percentage of sample chosen when doing Monte Carlo experiments
maxLoop = 1000; % nLoop times of Monte Carlo experiments
ci = 0.95; % confidence interval
% 20 years of data from Jan. 2005 to Dec. 2024
timeseries = datetime(2005, 1, 1):calmonths(1):datetime(2024, 12, 1);
timeseries = timeseries';

% parallel computing setting    
CoreNum = 24; % set the number of cores used in parallel computing
if isempty(gcp('nocreate'))
    parpool(CoreNum);
end

f = waitbar(0, 'Loading...');
% read the Excel file of groundwater depth at wells
path_input = '.\GroundwaterDepth.xlsx';
raw = readmatrix(path_input);

% type of monitoring well
% if 10 < type <= 19, the well monitors unconfined aquifers
% type 11 represents the first unconfined aquifers
% type 12 represents the second unconfined aquifers, etc.
% if 20 < type <= 29, the well monitors confined aquifers
% type 21 represents the first confined aquifers
% type 22 represents the second confined aquifers, etc.
type = raw(:, 2);
% longitude and latitude of the monitoring well
lon = raw(:, 3);
lat = raw(:, 4);
% monthly depth to groundwater depth data
data = raw(:, 5:end);

% type of well selected
% when doing interpolation, we only include wells in the first layer of
% unconfined (type = 11) and confined (type = 21) aquifers
typeMask = type == 21;

% reproject to UTM zone 50N (EPSG: 32650)
projUTM = projcrs(32650);
[utm_x, utm_y] = projfwd(projUTM, lat, lon);

% load the mask for interpolation
% choose your region of interest here
path_map = '.\mask\ncp.tif'; % entire North China Plain
[map, R] = readgeoraster(path_map);
infoR = geotiffinfo(path_map);

map_x = linspace(R.XWorldLimits(1), R.XWorldLimits(2), size(map, 2));
map_y = linspace(R.YWorldLimits(2), R.YWorldLimits(1), size(map, 1));
map_res = R.CellExtentInWorldX;

% find the coorditate of pixels in the selected region
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

% ---------------------------------------------------------------------
% find the location of in-situ groundwater depth in the map
waitbar(0, f, 'Start interpolation...');

avgDph = zeros(size(timeseries, 1), 1);
% parfor is used for parallel computing
parfor step = 1:size(timeseries, 1)
    % find available wells in the given month
    wellMask = typeMask & ~isnan(data(:, step));
    well_oi = data(wellMask, step);
    wellX_oi = utm_x(wellMask);
    wellY_oi = utm_y(wellMask);
    
    % find the rows and columns of grids that availble wells
    % belong to in the interpolated map 
    wellLoc = [zeros(size(well_oi, 1), 2), well_oi];
    wellLoc(:, 2) = sum(map_x < wellX_oi + 0.5*map_res, 2); % column
    wellLoc(:, 1) = sum(map_y > wellY_oi - 0.5*map_res, 2); % row
    wellLoc = sortrows(wellLoc);
    
    % if at least two available wells fall into the same grid,
    % calculate the average depth of all wells within the grid
    wellLoc_merge = wellLoc*0;
    mergeMask = zeros(size(wellLoc, 1), 1);
    for i = 1:size(wellLoc, 1)-1
        if wellLoc(i, 1) == wellLoc(i+1, 1) && wellLoc(i, 2) == wellLoc(i+1, 2)
            mergeMask(i) = 1;
        end
    end
    tempDph = 0;
    tempCount = 0;
    count = 0;
    for i = 1:size(wellLoc, 1)
        if mergeMask(i) == 1
            tempCount = tempCount + 1;
            tempDph = tempDph + wellLoc(i, 3);
        elseif mergeMask(i) == 0 & tempCount > 0
            count = count + 1;
            tempCount = tempCount + 1;
            tempDph = tempDph + wellLoc(i, 3);
            wellLoc_merge(count, :) = [wellLoc(i, 1:2), tempDph/tempCount];
            tempDph = 0;
            tempCount = 0;
        else
            count = count + 1;
            wellLoc_merge(count, :) = wellLoc(i, :);
        end
    end
    wellLoc_merge(count+1:end, :) = [];
    
    % 2D nearest neighbour interpolation (equivalent of Thiessen polygons)
    intNN = griddatan(wellLoc_merge(:, 1:2), wellLoc_merge(:, 3),...
                        loc_map(:, 1:2), 'nearest');
    
    mapNN = double(map*0);
    count = 0;
    for i = 1:size(map, 1)
        for j = 1:size(map, 2)
            if map(i, j) == 0
                % for grid out of the region of interest, we label as -99
                mapNN(i, j) = -99;
            else
                count = count + 1;
                mapNN(i, j) = intNN(count);
            end
        end
    end
    
    % save the geotiff result of spatial interpolation using following code
    % charMonth = num2str(year(timeseries(step))*100+month(timeseries(step)));
    % path_out = ['.\dphMap_unconfined\', charMonth, '.tif'];
    % geotiffwrite(path_out, mapNN, R, 'GeoKeyDirectoryTag', infoR.GeoTIFFTags.GeoKeyDirectoryTag);
    
    % calculate the weighted average depth
    avgDph(step) = mean(mapNN(mapNN~=-99), 'all', 'omitmissing');
end

% Uncertainty analysis

% record results of Monte Carlo simulation
avgDph_loop = zeros(size(timeseries, 1), maxLoop);

% Except for randomly selecting nSample of the input wells,
% everything else is the same as the code before
for nLoop = 1:maxLoop
    waitbar(nLoop/maxLoop, f, ['nLoop = ', num2str(nLoop)]);
    parfor step = 1:size(timeseries, 1)
        wellMask = typeMask & ~isnan(data(:, step));
        well_oi = data(wellMask, step);
        wellX_oi = utm_x(wellMask);
        wellY_oi = utm_y(wellMask);
    
        % randomly select nSample of the input wells
        randomMask = rand(size(well_oi)) < nSample;
        well_oi = well_oi(randomMask);
        wellX_oi = wellX_oi(randomMask);
        wellY_oi = wellY_oi(randomMask);
        
        wellLoc = [zeros(size(well_oi, 1), 2), well_oi];
        wellLoc(:, 2) = sum(map_x < wellX_oi + 0.5*map_res, 2);
        wellLoc(:, 1) = sum(map_y > wellY_oi - 0.5*map_res, 2);
        wellLoc = sortrows(wellLoc);
        
        wellLoc_merge = wellLoc*0;
        mergeMask = zeros(size(wellLoc, 1), 1);
        for i = 1:size(wellLoc, 1)-1
            if wellLoc(i, 1) == wellLoc(i+1, 1) && wellLoc(i, 2) == wellLoc(i+1, 2)
                mergeMask(i) = 1;
            end
        end
        tempDph = 0;
        tempCount = 0;
        count = 0;
        for i = 1:size(wellLoc, 1)
            if mergeMask(i) == 1
                tempCount = tempCount + 1;
                tempDph = tempDph + wellLoc(i, 3);
            elseif mergeMask(i) == 0 & tempCount > 0
                count = count + 1;
                tempCount = tempCount + 1;
                tempDph = tempDph + wellLoc(i, 3);
                wellLoc_merge(count, :) = [wellLoc(i, 1:2), tempDph/tempCount];
                tempDph = 0;
                tempCount = 0;
            else
                count = count + 1;
                wellLoc_merge(count, :) = wellLoc(i, :);
            end
        end
        wellLoc_merge(count+1:end, :) = [];
        
        intNN = griddatan(wellLoc_merge(:, 1:2), wellLoc_merge(:, 3),...
                            loc_map(:, 1:2), 'nearest');
        
        mapNN = double(map*0);
        count = 0;
        for i = 1:size(map, 1)
            for j = 1:size(map, 2)
                if map(i, j) == 0
                    mapNN(i, j) = -99;
                else
                    count = count + 1;
                    mapNN(i, j) = intNN(count);
                end
            end
        end
    
        avgDph_loop(step, nLoop) = mean(mapNN(mapNN~=-99), 'all', 'omitmissing');
    end
end
close(f);

% calculate the upper and lower bound of confidence interval
sorted = sort(avgDph_loop, 2, 'ascend');
lowerLoop = sorted(:, floor(0.5*(1-ci)*maxLoop)+1 );
upperLoop = sorted(:, ceil(0.5*(1+ci)*maxLoop) );

% visualize the results
close all;
plot(timeseries, avgDph, 'r-', 'LineWidth', 1);
hold on;
fill([timeseries; flipud(timeseries)], [lowerLoop; flipud(upperLoop)],...
    'r', 'FaceAlpha', 0.3, 'EdgeColor','none');
legend('AvgDph', 'Uncertainty');
ylabel('Depth (m)');
set (gca,'YDir','reverse');
hold off;

% find the final result here:
% regionally average depth, lower bound of confidence interval, and
% upper bound of confidence interval
AttentionForResult = [avgDph, lowerLoop, upperLoop];

toc; % record the time used for the script
