% Author: Yuancheng Xu
% Affiliation: Prof. Di Long's group at Tsinghua University, Beijing, PRC
% Email: xuyc24@mails.tsinghua.edu.cn
% Last updated: 25/04/2025

% Spatial interpolation: Thiessen polygons (Nearest neighbour)
function mapNN = GW_Dph_NN(well_oi, wellX_oi, wellY_oi, map_x, map_y, map_res, map, loc_map)
    
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
    