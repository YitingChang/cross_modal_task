%% high-speed video analysis for whisking
% whisker position captured by DeepLabCut
% input: x, y label position in pixels and the likelihood for each frame
% per body part
% frame rate: 
% to-do-list:
% remove outliers in DLC
% 
% Hilbert transformation: extracting the complex components from a signal (such as amplitude and phase)
% Band-pass filter data first! Hilbert transform can only be applied
% narrow-band signal.


MouseName = 'YT081';
videoDir = ['E:\', MouseName, '\SessionExplore_', MouseName];

hsVideo_path = 'C:\DeepLabCut\test\11-30-39.000DeepCut_resnet50_crossmodalSep18shuffle1_1030000.csv';
tb = readtable(hsVideo_path);
interval = (60*7+55)*1000/(size(tb,1)-2); % ms

selected_bodypart = {'surMid_x' 'surMid_y' 'surEnd_x' 'surEnd_y'};
selected_bodypart_ind = [11 12 14 15];

outlier_ind = [];
for i = 1:length(selected_bodypart)
    bodypart(:,i) = cellfun(@(x) str2num(x), tb{3:end, selected_bodypart_ind(i)}, 'UniformOutput', false);

    % post-hoc: remove outliers (> 5std)
    bodypart_mean = mean(cell2mat(bodypart(:,i)));
    bodypart_std = std(cell2mat(bodypart(:,i)));
    new_outlier_ind = cellfun(@(x) abs(x-bodypart_mean)/bodypart_std >5, bodypart(:,i), 'UniformOutput', false);
    outlier_ind = [outlier_ind new_outlier_ind]; 
end

% whisking threshold
k = theta>2;
whisking = find(k);
whisking_tp = (whisking/92465*475)/60; 

a = sum(cell2mat(outlier_ind), 2); 
b = cellfun(@(x) x == 0, num2cell(a), 'UniformOutput', false); % get outlier index
ind = repmat(cell2mat(b),1,4);
bodypart_02 = reshape(bodypart(ind), [length(bodypart(ind))/4,4]);


delta_x = cell2mat(bodypart_02(:,3)) - cell2mat(bodypart_02(:,1));
delta_y = cell2mat(bodypart_02(:,4)) - cell2mat(bodypart_02(:,2)); 
[theta, rho] = cart2pol(delta_x, delta_y);
delta_x = cell2mat(bodypart(:,3)) - cell2mat(bodypart(:,1));
delta_y = cell2mat(bodypart(:,4)) - cell2mat(bodypart(:,2)); 
[theta, rho] = cart2pol(delta_x, delta_y);
plot(theta);

% Butterworth bandpass (8-30Hz)
sampleRate = 1000/interval; % Hz
[b,a]=butter(2,[8 30]/sampleRate/2,'bandpass'); % digital or analog filters?
theta_02 = filtfilt(b, a, theta);

h = hilbert(theta_02);

abs(h)
phase = angle(h);

plot(phase)

subplot(2,1,1)
scatter(cell2mat(bodypart_02(:, 1)), cell2mat(bodypart_02(:, 2))); 
xlim([260 380])
subplot(2,1,2)
scatter(cell2mat(bodypart_02(:, 3)), cell2mat(bodypart_02(:, 4))); 
xlim([260 380])



% Get phase and frequency
% example:
t = 0:0.01:100;
y = sin(t);
hilbert(y);      
