%% Complete data array for each recSite
% Concatenate sessions 
% Fill with nan
% Get firing rate average

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth_lickOnset\data_array.mat')
% load('E:\Data_array\AllTrials_100bin_noSmooth\data_array.mat')

% Import setting    
B = 2; % number of contexts: respond-to-touch and respond-to-light blocks
S = 2; % number of stimuli: tactile and visual stimuli
D = 3; % number of decisions: lick (including right and left lick) and no lick

% recSites = {'left S1', 'left S2', 'left wM2', 'left AMM', 'left ALM', 'right ALM'};
recSites = {'left S1', 'left S2', 'left wM2', 'left ALM'};
% recSites = {'left AMM'};

for i = 1: length(recSites)
    u = data_all(strcmp(data_all.recSite, recSites{i}), :);
    maxTrialNum = max(cell2mat(u{:,5}),[],'all');
    Xtrial = [];
    trialNum = [];
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1};
        trialNum_session = u{session,5}{1,1};
        N = size(firingRates_session,1); % number of neurons
        for n = 1:N
            for b = 1:B
                for s = 1:S
                    for d = 1:D
                        firingRates_session(n,b,s,d,:,trialNum_session(n,b,s,d)+1:maxTrialNum) = nan;
                    end
                end
            end
        end         
        Xtrial = [Xtrial;firingRates_session];
        trialNum = [trialNum; trialNum_session];
    end
    firingRates{i,1} = recSites{i};
    firingRates{i,2} = Xtrial;
    firingRates{i,3} = trialNum;
end


%% grand average
% zScore:(x- mean(baseline))/std(baseline) 
%
% normalization: (x- mean(baseline))/std(baseline) vs. (x -mean(x))/range(x)
% (x- mean(baseline))/std(baseline) - (1) baseline firing rate is around zero (2) reflect actual firing rates (ex: low firing rate
% neurons remain low) - better for grand average    
%
% (x -mean(x))/range(x)- (1) normalized firing rate <= 1 but baseline firing rate may be negative (2) equal weights on all
% neurons regardless of their firing rates (ex: low firing rate neurons are scaled up) (3) smoother trace 
% - better for PCA because PCA is sensitivie to outliers
% 

mainDir = 'E:\CM_NeuralActivity_Analysis\zScore';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
baselineWindow = [-1 0];
baseline_startBin = (baselineWindow(1)+1)*(1/bin) + 1;
baseline_endBin = (baselineWindow(2)+1)*(1/bin);
stimWindow = [0 0.15];
nboot = 1000;
plotWindow = [-1 2.5];

tts = {{1 1 1} {1 1 3} {2 1 1:2} {2 1 3};...
    {2 2 2} {2 2 3} {1 2 1:2} {1 2 3};}; 

tt_names = {'Hit' 'Miss' 'FA' 'CR'};
tt_colors= {[0 0 1] [0.15 0.15 0.15] [0.4 0.8 0.1] [1 0 0]};
titles = {'Tactile stimulus' 'Visual stimulus'};


for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
    Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
    Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    
    figure('Position',[0 0 800 300])
    for s=1:S
        subplot(1,2,s)
        for tt=1:length(tts)
            Xtrial_tt = Xtrial(:,tts{s,tt}{1},tts{s,tt}{2},tts{s,tt}{3},:,:);
            % zScore
            X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials and different types of FA
            X_tt_normalized = (X_tt - Baseline_mean)./Baseline_std; 
            % bootstrapping for each time bin 
            X_tt_normalized = num2cell(X_tt_normalized,1); 
            bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), X_tt_normalized, 'UniformOutput', false);      
            bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
            bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)]', bootstat_sorted,'UniformOutput', false);
            bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           

            plot(time,bootstrap_summary(1,:),'Color', tt_colors{tt}); hold on;
            fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
                tt_colors{tt},'FaceAlpha', 0.3, 'Linestyle', 'none');
            plot(stimWindow, [1.2 1.2],'k','Linewidth',2)
            set(gca, 'box','off','TickDir','out')
            xlim(plotWindow);ylim([-0.2 1.2]);
            ylabel('Normalized firing rates')
            xlabel('Time from stimulus onset (s)')
            title(titles{s})
        end   
    end
    sgtitle(recSites{i}) 
    % Save figure  
    zScoreFigPath = fullfile(mainDir,'Figures_grandAverage',[recSites{i} '_zScore_all_trials_-1to2p5']);
    print(zScoreFigPath,'-dpdf','-painters','-loose');
end

%% Heatmap
mainDir = 'E:\CM_NeuralActivity_Analysis\zScore';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
baselineWindow = [-1 0];
baseline_startBin = (baselineWindow(1)+1)*(1/bin) + 1;
baseline_endBin = (baselineWindow(2)+1)*(1/bin);
stimWindow = [0 0.15];
plotWindow = [-0.25 1];
sortWindow = [0 0.5];
sort_startBin = (sortWindow(1)+1)*(1/bin) + 1;
sort_endBin = (sortWindow(2)+1)*(1/bin);

tts = {{1 1 1} {1 1 3} {1 2 1:2} {1 2 3};...
    {2 1 1:2} {2 1 3} {2 2 2} {2 2 3}}; 
tt_names = {'tHit' 'tMiss' 'vFA' 'vCR';...
    'tFA' 'tCR' 'vHit' 'vMiss'};
titles = {'Respond-to-touch block' 'Respond-to-light block'};


for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
    Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
    Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    
    figure('Position',[0 0 800 1000])
    for b=1:B            
        for tt=1:length(tts)
            Xtrial_tt = Xtrial(:,tts{b,tt}{1},tts{b,tt}{2},tts{b,tt}{3},:,:);
            % zScore
            X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials and different types of FA 
            X_tt_normalized = (X_tt - Baseline_mean)./Baseline_std; 
            %sorted by tHit 0~0.5 s
            if b==1 && tt==1
                [~,sort_idx] = sortrows(nanmean(X_tt_normalized(:,sort_startBin:sort_endBin),2),'descend');
            end
            X_normalized_sorted = X_tt_normalized(sort_idx,:);
            % plot
            subplot(2,4,(b-1)*4+tt)
            xTickLables = [repmat(' ',99,1); '0'; repmat(' ',99,1);'1'; repmat(' ',99,1); '2'; repmat(' ',50,1)]; 
            yTickLables = [{num2str(N)}; repmat({' '},N-2,1); {'1'}]; 
            h = heatmap(X_normalized_sorted,'Colormap',jet, 'ColorLimits',[-2 3],'Xlimits', [1,length(time)],...
                'ColorbarVisible', 'off','GridVisible','off',...
                'Title',tt_names{b,tt},'XDisplayLabels',xTickLables);
            if tt==1
                h.YLabel = 'Unit';
                h.YDisplayLabels= yTickLables;
            else
                h.YDisplayLabels= repmat({' '},N,1); 
            end
            if b==2
                h.XLabel = 'Time (s)';
            end
            if tt ==4
                h.ColorbarVisible = 'on';
            end
        end   
    end
    sgtitle(recSites{i}) 
    % Save figure  
    zScoreFigPath = fullfile(mainDir,'Figures_heatMap',[recSites{i} '_heatMap']);
    print(zScoreFigPath,'-dpdf','-painters','-loose');
end

%% tHit vs. tCR
% zScore (0~0.25)
mainDir = 'E:\CM_NeuralActivity_Analysis\zScore';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
baselineWindow = [-1 0];
baseline_startBin = (baselineWindow(1)+1)*(1/bin) + 1;
baseline_endBin = (baselineWindow(2)+1)*(1/bin);
analysisWindow = [0 0.25];
analysis_startBin = (analysisWindow(1)+1)*(1/bin) + 1;
analysis_endBin = (analysisWindow(2)+1)*(1/bin);
recSites = {'left S1', 'left S2', 'left wM2', 'left ALM'};
tts = {{1 1 1} {1 3 2}}; 

figure('Position',[0 0 800 800])  
for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
    Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
    Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    
    for tt=1:length(tts)
        Xtrial_tt = Xtrial(:,tts{tt}{1},tts{tt}{2},tts{tt}{3},analysis_startBin:analysis_endBin,:);
        % zScore
        X_tt = nanmean(Xtrial_tt(:,:),2); % average across trial types and trials 
        X_tt_normalized(:,tt) = (X_tt - Baseline_mean)./Baseline_std;    
    end
    subplot(2,2,i)
    scatter(X_tt_normalized(:,1),X_tt_normalized(:,2)); hold on;
    line([-10 15],[-10 15],'Color','black','LineStyle','--')
    xlim([-1 5]);ylim([-1 5]);
    xlabel('tHit z-score')
    ylabel('tCR z-score')                
    title(recSites{i}) 
    clearvars X_normalized
end

% Save figure  
tCorrectFigPath = fullfile(mainDir,'Figures_tCorrect','zScore_tCorrect');
print(tCorrectFigPath,'-dpdf','-painters','-loose');

%% Grand average for PCA version
% (x -mean(x))/range(x)

% timeWindow for analysis
bin = 0.01;
timeWindow = [-0.25, 0.25];
time = timeWindow(1)+bin:bin:timeWindow(2);
pre_timeBin = (timeWindow(1)+1)*(1/bin)+1;
post_timeBin = (timeWindow(2)+1)*(1/bin);
% T = post_timeBin - pre_timeBin +1; % number of time points
stimWindow = [0 0.15];
stimOn_Bin = (stimWindow(1)-timeWindow(1))*(1/bin);
stimOff_Bin = (stimWindow(2)-timeWindow(1))*(1/bin)-1;
componentsToPlot = 3;
% componentsToAnalyze = 5;

% select trial types for PCA
stimulus = [1 1 2 2]; 
decision = [1 3 2 3];
context = [1 2 2 1];
colors = {[0 0 1] [0 1 1] [1 0 0] [1 0 1]};  % rows: tactile stimulus, visual stimulus; columns: touch block, light block  
LineStyles = {'-x' '-x' '-x' '-x'};

for i = 1: length(recSites)
    X_tt = firingRates{i,3};
    for tt = 1:length(stimulus)
        firingRatesAverage_tt = squeeze(X_tt(:,stimulus(tt),decision(tt),context(tt),pre_timeBin:post_timeBin));
        firingRatesAverage(:,:,tt) = firingRatesAverage_tt;
    end
    X_mean = nanmean(firingRatesAverage(:,:),2);
    X_range = range(firingRatesAverage(:,:),2);
    firingRatesAverage_normalized = (firingRatesAverage - X_mean)./X_range;
    % Grand average figure
    figure('Position',[0 0 600 600])
    for tt = 1:length(stimulus)
        firingRatesAverage_normalized_tt = firingRatesAverage_normalized(:,:,tt);
        grandAverage = mean(firingRatesAverage_normalized_tt);
        plot(time,grandAverage,LineStyles{tt},'Color',colors{tt}, 'LineWidth', 1);hold on;
    end
    xline(time(stimOn_Bin));

    clearvars firingRatesAverage 
end

%% Grand average
% zScore:(x- mean(baseline))/std(baseline) 
% correct trials
% Chang et al. Figure 1g

% normalization: (x- mean(baseline))/std(baseline) vs. (x -mean(x))/range(x)
% (x- mean(baseline))/std(baseline) - (1) baseline firing rate is around zero (2) reflect actual firing rates (ex: low firing rate
% neurons remain low) - better for grand average    
%
% (x -mean(x))/range(x)- (1) normalized firing rate <= 1 but baseline firing rate may be negative (2) equal weights on all
% neurons regardless of their firing rates (ex: low firing rate neurons are scaled up) (3) smoother trace 
% - better for PCA because PCA is sensitivie to outliers
% 

mainDir = 'E:\CM_NeuralActivity_Analysis\zScore';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
baselineWindow = [-1 0];
baseline_startBin = (baselineWindow(1)+1)*(1/bin) + 1;
baseline_endBin = (baselineWindow(2)+1)*(1/bin);
stimWindow = [0 0.15];
nboot = 1000;
plotWindow = [-0.25 1];

tts = {{1 1 1} {2 1 3};...
    {2 2 2} {1 2 3};}; 
tt_names = {'Hit' 'CR'};
tt_colors = {[0 0 1] [0 1 1]; [1 0 0] [1 0 1]};
titles = {'Tactile stimulus' 'Visual stimulus'};


for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
    Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
    Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    
    figure('Position',[0 0 800 400])
    for s=1:S
        subplot(1,2,s)
        for tt=1:length(tts)
            Xtrial_tt = Xtrial(:,tts{s,tt}{1},tts{s,tt}{2},tts{s,tt}{3},:,:);
            % zScore
            X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials and different types of FA  
            X_tt_normalized = (X_tt - Baseline_mean)./Baseline_std; 
            % bootstrapping for each time bin 
            X_tt_normalized = num2cell(X_tt_normalized,1); 
            bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), X_tt_normalized, 'UniformOutput', false);      
            bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
            bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)]', bootstat_sorted,'UniformOutput', false);
            bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           

            plot(time,bootstrap_summary(1,:),'Color', tt_colors{s,tt}); hold on;
            fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
                tt_colors{s,tt},'FaceAlpha', 0.3, 'Linestyle', 'none');
            plot(stimWindow, [1.2 1.2],'k','Linewidth',2)
            set(gca, 'box','off','TickDir','out')
            xlim(plotWindow);ylim([-0.2 1.2]);
            ylabel('Normalized firing rates')
            xlabel('Time from stimulus onset (s)')
            title(titles{s})
        end   
    end
    sgtitle(recSites{i}) 
    % Save figure  
    zScoreFigPath = fullfile(mainDir,'Figures_grandAverage',[recSites{i} '_zScore_correct_trials']);
    print(zScoreFigPath,'-dpdf','-painters','-loose');
end


%% Heatmap (tBlock vs. vBlock)
% initial activity of correct trials 
% zScore:(x- mean(baseline))/std(baseline) 
% initial activity: -100 ~ 0 ms
% Chang et al., Fig 3b

% Import setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\zScore\Figures_heatmap_initial_activity';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
B=2;
bin = 0.1;
timeWindow = [-1 2.5];
initialBin = (0-timeWindow(1))*(1/bin);
tts = {{1 1 1} {1 2 3};... % tHit, vCR
    {2 1 3} {2 2 2}};  % tCR, vHit
titles = {'Respond-to-touch block' 'Respond-to-light block'};
% Info for example unit
example_units = {'left wM2', 19; 'left ALM', 108};

figure('Position',[0 0 600 400]) 
for i = 1: length(recSites)
    recSite = recSites{i};
    Xtrial = firingRates{strcmp(firingRates(:,1),recSite),2};
    N = size(Xtrial,1);
    Xtrial_initial = squeeze(Xtrial(:,:,:,:,initialBin,:));
    X_initial_mean = nanmean(Xtrial_initial(:,:),2);
    X_initial_std = nanstd(Xtrial_initial(:,:),0,2);
    Ztrial_initial = (Xtrial_initial - X_initial_mean)./(X_initial_std+0.01);
    
    for b=1:B  
        Ztrial_block = [];
        for tt=1:length(tts)
            Ztrial_tt = squeeze(Ztrial_initial(:,tts{b,tt}{1},tts{b,tt}{2},tts{b,tt}{3},:));
            Ztrial_block = [Ztrial_block Ztrial_tt];
        end
        Z_initial(1:N,b) = nanmean(Ztrial_block,2); % average across trial types and trials 
    end   
    [Z_initial_sorted, index] = sortrows(Z_initial,1,'descend'); % sorted by initial activity during respond-to-touch blocks

    % plot     
    subplot(1,4,i)   
    xTickLables = {'touch' 'light'};
    isExample = cell2mat(cellfun(@(x) strcmp(x,recSite),example_units(:,1), 'UniformOutput', false));
    if sum(isExample)>0
        % Label example units in the heatmap of initial activity
        unit_number = example_units{isExample,2};
        isUnit = index == unit_number;
        Unit_idx = find(isUnit==1);
        yTickLables = [{num2str(N)}; repmat({' '},Unit_idx-2,1); {'>'}; repmat({' '},N-Unit_idx-1,1); {'1'}]; 
    else
        yTickLables = [{num2str(N)}; repmat({' '},N-2,1); {'1'}]; 
    end
    h = heatmap(Z_initial_sorted,'Colormap',jet,'ColorLimits',[-0.6 0.6],...
        'ColorbarVisible', 'off','GridVisible','off',...
        'Title', recSites{i},'XDisplayLabels',xTickLables,'YDisplayLabels',yTickLables);
    if i==1
        h.YLabel = 'Unit';
    end
    if i==4
       h.ColorbarVisible = 'on';
    end  
    heatmap_idx{i,1} = recSite;
    heatmap_idx{i,2} = index;
    clear Xtrial Z_initial 
end           
sgtitle('Normalized initial activity (-0.1-0 s)')
% Save figure  
zScoreFigPath = fullfile(mainDir,'Normalized_initial_activity_-p1-0_labeled');
print(zScoreFigPath,'-dpdf','-painters','-loose');
% Save index of heat map
idxPaths = fullfile(mainDir,'Index of heat map');
save(idxPaths, 'heatmap_idx');

%% grand average for different FA types
% zScore:(x- mean(baseline))/std(baseline) 

mainDir = 'E:\CM_NeuralActivity_Analysis\zScore';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
baselineWindow = [-1 0];
baseline_startBin = (baselineWindow(1)+1)*(1/bin) + 1;
baseline_endBin = (baselineWindow(2)+1)*(1/bin);
stimWindow = [0 0.15];
nboot = 1000;
plotWindow = [-1 2.5];

tts = {{1 1 1} {2 1 1} {2 1 2};...
    {2 2 2} {1 2 2} {1 2 1}}; 

tt_colors= {[0 0 1] [0 0.6 1] [1 0.8 0]; 
            [1 0 0] [1 0.6 0] [0 0.8 1]};
titles = {'Tactile stimulus' 'Visual stimulus'};


for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
    Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
    Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    
    figure('Position',[0 0 800 300])
    for s=1:S
        subplot(1,2,s)
        for tt=1:length(tts)
            Xtrial_tt = Xtrial(:,tts{s,tt}{1},tts{s,tt}{2},tts{s,tt}{3},:,:);
            % zScore
            X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials
            X_tt_normalized = (X_tt - Baseline_mean)./Baseline_std; 
            % bootstrapping for each time bin 
            X_tt_normalized = num2cell(X_tt_normalized,1); 
            bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), X_tt_normalized, 'UniformOutput', false);      
            bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
            bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)]', bootstat_sorted,'UniformOutput', false);
            bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           

            plot(time,bootstrap_summary(1,:),'Color', tt_colors{s,tt}); hold on;
            fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
                tt_colors{s,tt},'FaceAlpha', 0.3, 'Linestyle', 'none');
            plot(stimWindow, [1.2 1.2],'k','Linewidth',2)
            set(gca, 'box','off','TickDir','out')
            xlim(plotWindow);ylim([-0.2 1.2]);
            ylabel('Normalized firing rates')
            xlabel('Time from stimulus onset (s)')
            title(titles{s})
        end   
    end
    sgtitle(recSites{i}) 
    % Save figure  
    zScoreFigPath = fullfile(mainDir,'Figures_grandAverage',[recSites{i} '_zScore_FA_types_-1to2p5']);
    print(zScoreFigPath,'-dpdf','-painters','-loose');
end

%% grand average for all tactile trial types 
% zScore:(x- mean(baseline))/std(baseline) 


mainDir = 'E:\CM_NeuralActivity_Analysis\zScore';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
baselineWindow = [-1 0];
baseline_startBin = (baselineWindow(1)+1)*(1/bin) + 1;
baseline_endBin = (baselineWindow(2)+1)*(1/bin);
stimWindow = [0 0.15];
nboot = 1000;
plotWindow = [-1 2.5];

% tts = {{1 1 1} {2 1 1} {2 1 2} {};...
%     {1 1 3} {2 1 3} {[]}}; 

tt_colors= {[0 0 1] [1 0.4 0]; [0.5 0.2 0.5] [1 0.8 0]; [0.5 0.5 0.5] [0 1 1]};
% titles = {'Tactile + lick' 'Tactile + no lick'};

figure('Position',[0 0 900 500])
for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
    Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
    Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    for d=1:D
        subplot(3,4,(d-1)*4+i)
        for b=1:B
            Xtrial_tt = Xtrial(:,b,1,d,:,:);
            % zScore
            X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials
            X_tt_normalized = (X_tt - Baseline_mean)./Baseline_std; 
            % remove units that do not have certain trial types 
            isUnit = sum(isnan(X_tt_normalized),2) == 0;
            X_tt_normalized = num2cell(X_tt_normalized(isUnit,:),1); 
            % bootstrapping for each time bin
            bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), X_tt_normalized, 'UniformOutput', false);      
            bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
            bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)]', bootstat_sorted,'UniformOutput', false);
            bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           
            plot(time,bootstrap_summary(1,:),'Color', tt_colors{d,b}); hold on;
            fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
                tt_colors{d,b},'FaceAlpha', 0.3, 'Linestyle', 'none');
            plot(stimWindow, [1.2 1.2],'k','Linewidth',2)
        end
        set(gca, 'box','off','TickDir','out')
        xlim(plotWindow);ylim([-0.2 1.2]);
        ylabel('Normalized firing rates')
        xlabel('Time from stimulus onset (s)')
    end
end
% Save figure  
zScoreFigPath = fullfile(mainDir,'Figures_grandAverage','zScore_tactile_trials_-1to2p5');
print(zScoreFigPath,'-dpdf','-painters','-loose');


%% Permutation test 
% O'Connor et al., Neuron, 2010 (Fig. S7D)
% Determine whether the mean PSTHs for early and lare trials for each neuron were siginificant, as measured by the Euclidean
% distance between them (Foffani and Moxon, 2004; Sandler, 2008)
% clearvars -except firingRates numOfTrials trialNum_tBlock trialNum_vBlock transition
bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow = [0 0.25];
binStart =(analysisWindow(1)-timeWindow(1))/bin+1;
binEnd = (analysisWindow(2)-timeWindow(1))/bin;
N = size(firingRates,1); % number of neurons

% tts = {[1 1 1] [1 3 1] [1 1 2] [1 2 2] [1 3 2]}; % TTT (tHit), TTN (tMiss), VTT (tFA, rule error), VTV (tFA, compulsive), VTN (tCR)
% tt_names = {'tHit' 'tMiss' 'tFA (rule)' 'tFA (compulsive)' 'tCR'};

rng(17); % control random number generation
nShuffle = 1000; 
alpha = 0.05;
for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    trialNum = firingRates{strcmp(firingRates(:,1),recSites{i}),3};
    N = size(Xtrial,1);
%     Xtrial_baseline= Xtrial(:,:,:,:,baseline_startBin:baseline_endBin,:);
%     Baseline_mean = nanmean(Xtrial_baseline(:,:),2);
%     Baseline_std = nanstd(Xtrial_baseline(:,:),0,2);
    for d=1:D
        for n=1:N
            trialNum_A = trialNum(n,1,1,d);
            trialNum_B = trialNum(n,2,1,d);
            if trialNum_A>0 & trialNum_B>0
                Xtrial_A = squeeze(Xtrial(n,1,1,d,binStart:binEnd,1:trialNum_A));
                Xtrial_B = squeeze(Xtrial(n,2,1,d,binStart:binEnd,1:trialNum_B));
                X_A = mean(Xtrial_A,2); % respond-to-touch block
                X_B = mean(Xtrial_B,2); % respond-to-light block
                distance = norm(X_A - X_B);
                Xtrial_AB = [Xtrial_A Xtrial_B];
                for shuffle=1:nShuffle
                    new_trialNum = randperm(trialNum_A+trialNum_B);
                    isA = new_trialNum <= trialNum_A;
                    Xtrial_A_shuffled = Xtrial_AB(:,isA);
                    Xtrial_B_shuffled = Xtrial_AB(:,~isA);
                    X_A_shuffled = mean(Xtrial_A_shuffled,2); 
                    X_B_shuffled = mean(Xtrial_B_shuffled,2); 
                    distance_shuffled(shuffle) = norm(X_A_shuffled - X_B_shuffled);
                end
                pValues(d,n) = sum(distance_shuffled>=distance)/length(distance_shuffled);% one tailed 
            else
                pValues(d,n) = nan;
            end 
        end
        N_adjusted = sum(~isnan(pValues(d,:)));
        alpha_adjusted = alpha/N_adjusted; % Bonferroni correction
        isSig = pValues(d,:) < alpha_adjusted;
        summary_recSite(d,1) = N;
        summary_recSite(d,2) = N_adjusted;
        summary_recSite(d,3) = sum(isSig);
        summary_recSite(d,4) = sum(isSig)/N_adjusted;
    end
    summary{i,1} = recSites{i};
    summary{i,2} = summary_recSite;
    clear Xtrial trialNum pValues summary_recSite
end

