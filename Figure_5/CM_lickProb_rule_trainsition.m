%% Plot lick direction probability during rule trainsition
% Chang et al Fig 5a

% Get right and left lick probability during rule transition across mice
load('E:\Cross-Modal_Project\Revision\rule_transition\data_array.mat')
mainDir = 'E:\Cross-Modal_Project\Revision\rule_transition';

recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
u = data_all(ismember(data_all.recSite, recSites), :);
MouseNames = unique(u.MouseName);
B=2; % respond-to-touch and respond-to-light blocks
L=2; % right and left licks
lick_colors = {'b' 'r'};
block_titles = {'Light to Touch' 'Touch to Light'};
% Get 10 trials before block switch and 15 trials after block switch
transWindow = [10 15];
numTrans = transWindow(2) + transWindow(1) + 1; % number of trials points
for m=1:length(MouseNames)
    MouseName = MouseNames(m);
    v = u(strcmp(u.MouseName,MouseName),:); 
    for session = 1:size(v,1)
        lick(1,1,1:numTrans,session) = mean(v{session,4}{1,1},2); % trans_rLick_tBlock
        lick(1,2,1:numTrans,session) = mean(v{session,5}{1,1},2); % trans_lLick_tBlock
        lick(2,1,1:numTrans,session) = mean(v{session,6}{1,1},2); % trans_rLick_vBlock
        lick(2,2,1:numTrans,session) = mean(v{session,7}{1,1},2); % trans_lLick_vBlock
    end
    lick_all(1:2,1:L,1:numTrans,m) = mean(lick,4);
end

figure('Position',[0 0 800 400])
for b=1:B
    subplot(1,2,b)
    for l=1:L
        lick_tt = squeeze(lick_all(b,l,:,:));
        x = -transWindow(1):transWindow(2);
        y = mean(lick_tt,2);
        plot(x,y,lick_colors{l}); hold on;
    end
    title(block_titles{b})
    set(gca, 'box','off', 'TickDir','out')
    ylabel('Probability')
    xlabel('Trial number from block switch')
    xlim([-5 15])
    ylim([0 1])
end
% Save figure
FigPath = fullfile(mainDir,'lick_prob');
print(FigPath,'-dpdf','-painters','-loose');
        
%% Get right and left lick probability during rule transition across mice (Bootstrapping)
load('E:\Cross-Modal_Project\Revision\rule_transition\data_array.mat')
mainDir = 'E:\Cross-Modal_Project\Revision\rule_transition';

recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
u = data_all(ismember(data_all.recSite, recSites), :);
MouseNames = unique(u.MouseName);
B=2; % respond-to-touch and respond-to-light blocks
L=2; % right and left licks
lick_colors = {'b' 'r'};
block_titles = {'Light to Touch' 'Touch to Light'};
% Get 10 trials before block switch and 15 trials after block switch
transWindow = [10 15];
numTrans = transWindow(2) + transWindow(1) + 1; % number of trials points
nboot = 1000;
rng(23)
for p = 1:nboot
    % resample mice
    isMouse = randi(length(MouseNames), length(MouseNames), 1);
    resampled_MouseNames = MouseNames(isMouse);
    for m=1:length(resampled_MouseNames)
        MouseName = resampled_MouseNames{m};
        v = u(strcmp(u.MouseName,MouseName),:); 
        % resample sessions
        isSession = randi(size(v,1),size(v,1),1);
        resampled_v = v(isSession,:);
        for session=1:size(resampled_v,1)
            % resample transitions
            num_tBlock = size(resampled_v{session,4}{1,1},2);
            num_vBlock = size(resampled_v{session,6}{1,1},2);
            istBlock = randi(num_tBlock,num_tBlock,1);
            isvBlock = randi(num_vBlock,num_vBlock,1);
            
            lick(1,1,1:numTrans,session) = mean(resampled_v{session,4}{1,1}(:,istBlock),2); % trans_rLick_tBlock
            lick(1,2,1:numTrans,session) = mean(resampled_v{session,5}{1,1}(:,istBlock),2); % trans_lLick_tBlock
            lick(2,1,1:numTrans,session) = mean(resampled_v{session,6}{1,1}(:,isvBlock),2); % trans_rLick_vBlock
            lick(2,2,1:numTrans,session) = mean(resampled_v{session,7}{1,1}(:,isvBlock),2); % trans_lLick_vBlock
            clearvars istBlock isvBlock
        end
        lick_all(1:B,1:L,1:numTrans,m) = mean(lick,4);
        clearvars lick resampled_v
    end
    lick_bootstrapped(1:B,1:L,1:numTrans,p) = mean(lick_all,4);
    clearvars lick_all resampled_MouseNames
end

% Save bootstrapping results 
dataSavePath = fullfile(mainDir, 'lick_bootstrapped');
save(dataSavePath, 'lick_bootstrapped','-v7.3');

% Plot
x = -transWindow(1):transWindow(2);
x_bw = [x, fliplr(x)];
figure('Position',[0 0 800 400])
for b=1:B
    subplot(1,2,b)
    for l=1:L
        lick_tt = num2cell(squeeze(lick_bootstrapped(b,l,:,:)),2);
        lick_sorted = cellfun(@(x) sort(x), lick_tt, 'UniformOutput', false);
        bootstrap_summary = cellfun(@(x) [mean(x), x(round(0.025*nboot)), x(round(0.975*nboot))], lick_sorted,'UniformOutput', false);
        bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI  
        plot(x,bootstrap_summary(:,1),'Color', lick_colors{l}); hold on;
        fill(x_bw, [bootstrap_summary(:,2)', fliplr(bootstrap_summary(:,3)')],...
            lick_colors{l},'FaceAlpha', 0.3, 'Linestyle', 'none'); hold on;
    end
    title(block_titles{b})
    set(gca, 'box','off', 'TickDir','out')
    ylabel('Probability')
    xlabel('Trial number from block switch')
    xlim([-5 15])
    ylim([0 1])
end
                 
% Save figure
FigPath = fullfile(mainDir,'lick_prob_bootstrapped');
print(FigPath,'-dpdf','-painters','-loose');        
