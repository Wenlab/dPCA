%% load data and set parameters
path_working = 'D:\\Users\\Ali_Li\\Desktop\\to_daguang_0619\\';
load(fullfile(path_working,'ca_trace.mat'));
fps = 50;
T_beat_min = fps/1000.0*29;
T_beat_max = fps/1000.0*52;
T_bout_max = fps/1000.0*230;
curve_centerlines_interpolated_smooth = curvature;
curve_centerlines_detrend = curve_centerlines_interpolated_smooth;
length_extent = 0;
num_curvs_res = 8;

fps_fluore = 1;

%% behavior analysis: detect large and small bouts

%% Bending energy
bending_energy = sum(curve_centerlines_detrend.*curve_centerlines_detrend, 2);
mean_bending_energy = mean(bending_energy,'omitnan');
std_bending_energy = std(bending_energy,'omitnan');
% figure; plot(bending_energy);
% hold on;
% plot(2*std_bending_energy*ones(length(bending_energy),1) + mean_bending_energy);
% hold off;
% title('bending energy');

%% Kinetic energy
% abs_angle = cumsum(curve_centerlines_interpolated_smooth,2);
% abs_disp = cumsum(abs_angle,2);
% we should smooth the displacement first, then calculate velocity
% velocity = diff(abs_disp, 1, 1);
% kinetic_energy = sum(velocity.*velocity, 2);

angular_velocity = diff(curve_centerlines_detrend, 1, 1);
kinetic_energy = sum(angular_velocity.*angular_velocity, 2);
mean_kinetic_energy = mean(kinetic_energy,'omitnan');
std_kinetic_energy = std(kinetic_energy,'omitnan');

% figure; plot(kinetic_energy);
% hold on;
% plot(std_kinetic_energy*ones(length(kinetic_energy),1) + mean_kinetic_energy);
% hold off;
% title('kinetic energy');

%% detect bouts
normalized_bending_energy = (bending_energy-mean_bending_energy)/std_bending_energy;
normalized_kinetic_energy = (kinetic_energy-mean_kinetic_energy)/std_kinetic_energy;
normalized_kinetic_energy = [0;normalized_kinetic_energy];
mechanical_energy = normalized_bending_energy + normalized_kinetic_energy;
mechanical_energy_smooth = smoothdata(mechanical_energy,'movmean',3);

mean_mechanical_energy = mean(mechanical_energy_smooth,'omitnan');
std_mechanical_energy = std(mechanical_energy_smooth,'omitnan');
% figure; plot(mechanical_energy_smooth);
% hold on;
% plot(normalized_bending_energy);
% plot(normalized_kinetic_energy);
% plot(std_mechanical_energy*ones(length(mechanical_energy_smooth),1) + median(mechanical_energy_smooth,'omitnan'));
% hold off;
% title('energy');

% mechanical_energy_smooth = smoothdata(mechanical_energy,'gaussian',3);
% figure; plot(mechanical_energy);
% hold on;
% plot(mechanical_energy_smooth);

GMModel = fitgmdist(mechanical_energy_smooth,2);
% flag_bout = cluster(GMModel,mechanical_energy_smooth);% 确定一下背景和bout的编号
% [M,I] = max(mechanical_energy_smooth);
% ID_bout = flag_bout(I);
% [M,I] = min(mechanical_energy_smooth);
% ID_interval = flag_bout(I);

% flag_bout_ori = flag_bout;
% flag_bout = zeros(length(flag_bout_ori),1);
% flag_bout(flag_bout_ori==ID_bout) = 1;

NumPoints = (max(mechanical_energy_smooth)-min(mechanical_energy_smooth))/(max(GMModel.mu)-min(GMModel.mu))*10; % 找能量的分布密度的第一个极小值点，作为bout的阈值
NumPoints2 = (max(mechanical_energy_smooth)-min(mechanical_energy_smooth))/(max(GMModel.mu)-min(GMModel.mu))*100;
[f,xi,bw] = ksdensity(mechanical_energy_smooth,'NumPoints',NumPoints);
[f2,xi2,bw2] = ksdensity(mechanical_energy_smooth,'NumPoints',NumPoints2);
% figure;
% plot(xi,f);
% title('distribution of energy');

% TF = islocalmin(f);
% for i=1:length(TF)
%     if (TF(i))
%         threshold_energy = xi(i);
%         break;
%     end
% end
% set energy threshold mannually
threshold_energy = 0.0144; % energy threshold for bout or interval
threshold_energy2 = 0.865; % energy threshold for large or small bout

flag_nan = isnan(mechanical_energy_smooth);
idx_n = find(~flag_nan);
idx_nan = find(flag_nan);
mechanical_energy_smooth(idx_nan) = interp1(idx_n, mechanical_energy_smooth(idx_n), idx_nan);
flag_bout = mechanical_energy_smooth > threshold_energy;
flag_bout_ori = flag_bout;

flag_bout = bwareaopen(flag_bout, 3);% 短于3的噪点会被清除
SE = strel('line',round(T_beat_min),90);% 距离小于T_beat_min的两个bouts会被合并
flag_bout = imclose(flag_bout,SE);
flag_bout = bwareaopen(flag_bout,round(T_beat_min));% 去除短于T_beat_min的连通区域，bout不短于T_beat_min帧（一个最短摆尾周期）。
SE2 = strel('line',round(T_beat_min + length_extent*2),90);% 将bout像两侧各扩展T_beat_min/2 + length_extent帧
flag_bout = imdilate(flag_bout,SE2);

% figure;
% plot(mechanical_energy_smooth,'-o','MarkerIndices',find(flag_bout_ori==1));
% title('flag bout original');

start_bout = [];
end_bout = [];
pre_flag_bout = false;
for i=1:length(flag_bout)
    present_flag_bout = flag_bout(i);
    if (~pre_flag_bout && present_flag_bout)
        start_bout = [start_bout, i];
    elseif (pre_flag_bout && ~present_flag_bout)
        end_bout = [end_bout, i-1];
    end
    pre_flag_bout = present_flag_bout;
end
if (size(end_bout,2)<size(start_bout,2)) % The last bout has not finished.
    start_bout(end) = [];
end

flag_bout_intact = true(length(start_bout),1);
for i=1:length(start_bout)
    flag_nan = isnan(curve_centerlines_interpolated_smooth(start_bout(i):end_bout(i),1:num_curvs_res));
    num_nan = sum(flag_nan,'all');
    if (num_nan>0)% curve_centerlines_interpolated_smooth中前num_curvs_res段有NaN的bout，去除
        flag_bout_intact(i) = false;
    end
end
start_bout(~flag_bout_intact) = [];
end_bout(~flag_bout_intact) = [];

% figure; 
% plot(mechanical_energy_smooth);
% hold on;
% scatter(start_bout,mechanical_energy_smooth(start_bout),23,'g');
% scatter(end_bout,mechanical_energy_smooth(end_bout),23,'r');
% title('bouts start end');

%% save bout detection results
num_bouts = length(start_bout);
bout_template = struct('start_bout', single(0), 'end_bout', single(0), 'duration', single(0), 'mid_bout_fluore', single(0), 'bout_type', single(0));
% bout_type: 0: interval; 1: small bout; 2: large bout.
bouts = repmat(bout_template,[num_bouts 1]);
for i=1:num_bouts
    bouts(i).start_bout = single(start_bout(i));
    bouts(i).end_bout = single(end_bout(i));
    bouts(i).duration = single(end_bout(i) - start_bout(i) + 1);
    mid_bout = (bouts(i).start_bout + bouts(i).end_bout)*0.5;
    bouts(i).mid_bout_fluore = single(round(mid_bout/fps*fps_fluore));
    if max(mechanical_energy_smooth(bouts(i).start_bout:bouts(i).end_bout))>threshold_energy2
        bouts(i).bout_type = single(2);
    else
        bouts(i).bout_type = single(1);
    end
end

save(fullfile(path_working, ['bout_detection_withHeading_extent', num2str(length_extent),'.mat']), 'bouts','flag_bout_ori','-v7.3');

%% show the bout detection result
% figure;
% plot(curve_centerlines_detrend);
% hold on; 
% for i=1:num_bouts
%     temp = ones(bouts(i).duration,1)*bouts(i).bout_type*0.1;
%     plot(bouts(i).start_bout:bouts(i).end_bout,temp,'LineWidth',3);
% end

