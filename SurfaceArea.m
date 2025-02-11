% Uploaded by YoungHye Judy Kwon (Dec 4, 2024) 

clear; clc;
% Add required paths
addpath('../cifti-matlab-master');

% Subject and network IDs
subID = [# # # # # #];
pmnID = [# # # # # #];
dnaID = [# # # # # #];
dnbID = [# # # # # #];
fpnaID = [# # # # # #];
fpnbID = [# # # # # #];
salID = [# # # # # #];
conID = [# # # # # #];

% Define regions and networks
region_lh = {'ant_lateral_lh', 'ant_medial_lh', 'post_lateral_lh', 'pos_medial_lh'};
region_rh = {'ant_lateral_rh', 'ant_medial_rh', 'post_lateral_rh', 'pos_medial_rh'};
network = {'dnaID', 'dnbID', 'pmnID', 'fpnaID', 'fpnbID', 'salID', 'conID'};
hemis = {'left', 'right'};
whole_vert_num = 327684;   % 327684 for fsaverage7

% Process each subject and network
for s = 1:numel(subID)
    for n = 1:numel(network)
        for i = 1:numel(region_lh) % Loop over regions
            % Load parcellation maps
            tpl = cifti_read(['/path/to/parcellation/map']);   % load dlabel file

            % Get current network IDs
            a=zeros(numel(tpl.cdata),1);
            cur_network=eval(network{n});
            a(find(tpl.cdata == cur_network(s)),1)=1;
            tpl.cdata = a;

            % divide anterior-posterior for SAL/PMN
            if n == 3 % if current network is SAL/PMN
                tmp_lh = cifti_read(['/path/to/vertex/file/']) % dlabel file for left hemisphere
                tmp_rh = cifti_read(['/path/to/vertex/file/']) % dlabel file for right hemisphere
                net_vert = [tmp_lh.cdata; tmp_rh.cdata]; % combine vertex IDs from both hemispheres
            else
                net_vert = a;
            end

            % Remove medial wall
            medialwall = cifti_read(['/path/to/medial/wall/']); % load medial wall indices (dlabel)
            mw_count = sum(medialwall.cdata);
            tmp = net_vert - medialwall.cdata;
            tmp(tmp == -1) = 0;

            % Get vertex indices for the region
            idx = find(tmp &  tpl.cdata(1:numel(tmp)));

            % Count vertices in left and right hemispheres
            count_vert{n}{s}(1, i) = length(find(idx < whole_vert_num / 2)); % Left hemisphere
            count_vert{n}{s}(2, i) = length(find(idx > whole_vert_num / 2)); % Right hemisphere
        end
    end
end

disp('Initial vertex counting done.');

% Compute percentages for each network
for s = 1:numel(subID)
    for n = 1:numel(network)
        if n == 3 % SAL/PMN (ant-pos division)
            num_vert{n}{s}(1, :) = sum(sum(count_vert{n}{s}(:, 1:2))) / (whole_vert_num - mw_count) * 100; % Anterior
            num_vert{n}{s}(2, :) = sum(sum(count_vert{n}{s}(:, 3:4))) / (whole_vert_num - mw_count) * 100; % Posterior
            num_vert{n}{s}(3, :) = sum(num_vert{n}{s}(1:2, :)); % Anterior + Posterior
        else % Other networks
            num_vert{n}{s}(1, :) = count_vert{n}{s}(1, 1) / (whole_vert_num - mw_count) * 100; % Left hemisphere
            num_vert{n}{s}(2, :) = count_vert{n}{s}(2, 1) / (whole_vert_num - mw_count) * 100; % Right hemisphere
            num_vert{n}{s}(3, :) = sum(num_vert{n}{s}(1:2, :)); % Both hemisphere
        end
    end
end

disp('Percentage calculation done.');

%% Statistical Analysis
disp('Performing statistical tests...');
% Unified SAL/PMN vs. other networks (plot1)
for n=1:numel(network)
    for s=1:numel(subID)
        % for plot1
        ydata(n,s) = num_vert{n}{s}(3);
    end
    ydata_mean(n,1) = mean(ydata(n,:))
end

% Perform ANOVA
[p,tbl,stats]=anova1(ydata');
results=multcompare(stats);

% Splited SAL/PMN vs. other networks (plot2)
for n=1:numel(network)
    for s=1:numel(subID)
        if n==3
            ydata{n}(s,1) = num_vert{n}{s}(1); %anterior
            ydata{n}(s,2) = num_vert{n}{s}(2); %posterior
        else
            ydata{n}(s) = num_vert{n}{s}(3);  %all
        end
    end
end

% Perform ANOVA
for y=1:numel(ydata)
    if y==3
        ydata2(y,:)=ydata{y}(:,1)';
        ydata2(y+1,:)=ydata{y}(:,2)';
    elseif y>3
        ydata2(y+1,:)=ydata{y};
    else
        ydata2(y,:)=ydata{y};
    end
end
ydata2=ydata2';
[p,tbl,stats]=anova1(ydata2);
results=multcompare(stats);

disp('Statistical analysis completed.');
