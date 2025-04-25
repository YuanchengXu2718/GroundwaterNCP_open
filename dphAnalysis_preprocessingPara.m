% Author: Yuancheng Xu
% Affiliation: Prof. Di Long's group at Tsinghua University, Beijing, PRC
% Email: xuyc24@mails.tsinghua.edu.cn
% Last updated: 15/01/2025

% Code for spatial interpolation and uncertainty analysis of regionally 
% average groundwater depth

tic; % record the time used for the script

% path_ref = "E:\2024groundwater_insitu\In_situ_writing_revise1\cumDeltaDph.xlsx";
% ref = readmatrix(path_ref, 'Range', 'Q2:Q241');

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

% type of well selected
% when doing interpolation, we only include wells in the first layer of
% unconfined (type = 11) and confined (type = 21) aquifers
typeMask = type == 11;

% reproject to UTM zone 50N (EPSG: 32650)
projUTM = projcrs(32650);
[utm_x, utm_y] = projfwd(projUTM, lat, lon);

% load the mask for interpolation
% choose your region of interest here
path_map = root + "\mask\ncp.tif"; % entire North China Plain
[map, R] = readgeoraster(path_map);
infoR = geotiffinfo(path_map);

map_x = linspace(R.XWorldLimits(1), R.XWorldLimits(2), size(map, 2));
map_y = linspace(R.YWorldLimits(2), R.YWorldLimits(1), size(map, 1));
map_res = R.CellExtentInWorldX;

intrplGAP_ = 1:5;
minLen_ = [12, 24, 36, 48, 60];
windowOL_ = [12, 24, 36, 48, 60];

res_para = zeros(size(dphdata, 2), length(intrplGAP_)*15);
input_para = zeros(3, length(intrplGAP_)*15);
% res_para_delta = zeros(size(dphdata, 2), length(interLen_));
col_res = 0;

dphdata_raw = dphdata;

% temporal interpolation (linear)
for intrplGAP_i = 1:length(intrplGAP_)
    intrplGAP = intrplGAP_(intrplGAP_i);
    dphdata = dphdata_raw;
    for i = 1:size(dphdata, 1)
        temInter = 0;
        for j = 1:size(dphdata, 2)
            if temInter > 0
                if ~isnan(dphdata(i, j)) && j == temInter + 1 % valid, valid
                    temInter = j;
                elseif ~isnan(dphdata(i, j)) && j <= temInter + intrplGAP + 1 % valid, nan (<= interLen), valid 
                    dphdata(i, temInter:j) = linspace(dphdata(i, temInter), dphdata(i, j), j - temInter + 1);
                    temInter = j;
                elseif ~isnan(dphdata(i, j)) && j > temInter + intrplGAP + 1 % valid, nan (> interLen), valid
                    temInter = j;
                end
            elseif ~isnan(dphdata(i, j))
                temInter = j;
            else
                continue
            end
        end
    end       
    
    dphdata_tempInter = dphdata;
    for minLen_i = 1:5
        for windowOL_i = 1:minLen_i
            dphdata = dphdata_tempInter;
            col_res = col_res + 1;
            waitbar(col_res/size(res_para, 2), f, "case = " + num2str(col_res));
            minLen = minLen_(minLen_i);
            outlier_window = windowOL_(windowOL_i);

            for i = 1:size(dphdata, 1)
                if sum(~isnan(dphdata(i, :))) < minLen
                    dphdata(i, :) = nan;
                end
                outlierMask = zeros(size(dphdata, 2), 1);
                for j = outlier_window:size(dphdata, 2)
                    dph_window = dphdata(i, j - outlier_window + 1 : j);
                    avg_window = mean(dph_window, 'omitmissing');
                    sigma_window = sqrt(var(dph_window, 'omitmissing'));
                    mask_window = dph_window > avg_window + 3*sigma_window |...
                                    dph_window < avg_window - 3*sigma_window;
                    outlierMask(j - outlier_window + 1 : j) = mask_window' + outlierMask(j - outlier_window + 1 : j);
                end
                outlierMask = outlierMask > 0.5 * outlier_window;
                dphdata(i, outlierMask) = nan;
            end
            
            for i = 2:size(data, 2)
                data(:, i) = dphdata(:, i) - dphdata(:, i-1);
            end   
            
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
            
            % find the location of in-situ groundwater depth in the map         
            deltaDph = zeros(size(timeseries, 1), 1);
            cumDph = zeros(size(timeseries, 1), 1);
            % parfor is used for parallel computing
            parfor step = 2:size(timeseries, 1)
                % find available wells in the given month
                wellMask = typeMask & ~isnan(data(:, step));
                well_oi = data(wellMask, step);
                wellX_oi = utm_x(wellMask);
                wellY_oi = utm_y(wellMask);
                
                mapNN = GW_Dph_NN(well_oi, wellX_oi, wellY_oi, map_x, map_y, map_res, map, loc_map);
                
                % calculate the weighted average delta depth
                deltaDph(step) = mean(mapNN(mapNN~=-99), 'all', 'omitmissing');
            end
            
            for step = 2:size(timeseries, 1)
                cumDph(step) = cumDph(step - 1) + deltaDph(step);
            end
            res_para(:, col_res) = cumDph;
            input_para(:, col_res) = [intrplGAP; minLen; outlier_window];
        end
    end
end

close(f);

sorted = sort(res_para, 2, 'ascend');
lowerLoop = sorted(:, 1);
upperLoop = sorted(:, end);

% visualize the results
close all;
figure(1);
% parameters used in the results: intrplGAP = 3, minLen = 48, windowOL = 24
ref = res_para(:, 38);
plot(timeseries, ref - mean(ref), '-r', 'LineWidth', 1, 'MarkerSize', 1);
hold on
fill([timeseries; flipud(timeseries)], [lowerLoop - mean(ref); flipud(upperLoop) - mean(ref)],...
    'r', 'FaceAlpha', 0.3, 'EdgeColor','none');
hold off
ylabel('Groundwater depth anomaly (m)');
set (gca,'YDir','reverse');

% find the final result here:
% regionally average depth, lower bound of confidence interval, and
% upper bound of confidence interval
AttentionForResult = [ref - mean(ref), lowerLoop - mean(ref), upperLoop - mean(ref)];
round(AttentionForResult, 2);

toc; % record the time used for the script
