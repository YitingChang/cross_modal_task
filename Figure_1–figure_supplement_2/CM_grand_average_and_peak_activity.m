% Revision 
% Get grand averages in firing rates instead of zScore
% Get distribution of peak activity to show why it is important to normalize neural activity before average across neurons

%% Complete data array for each recSite
% Concatenate sessions 
% Fill with nan
% Get firing rate average

% load data array
load('E:\Cross-Modal_Project\Data_array\AllTrials_10bin_Smooth_lickOnset\data_array.mat')
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



%% Grand average (firing rates in Hz) for correct trials 

mainDir = 'E:\Cross-Modal_Project\Revision\grand_average_and_peak_activity';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
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
    
    figure('Position',[0 0 800 400])
    for s=1:S
        subplot(1,2,s)
        for tt=1:length(tts)
            Xtrial_tt = Xtrial(:,tts{s,tt}{1},tts{s,tt}{2},tts{s,tt}{3},:,:);
            X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials and different types of FA  
            % bootstrapping for each time bin 
            X_tt = num2cell(X_tt,1); 
            bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), X_tt, 'UniformOutput', false);      
            bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
            bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)]', bootstat_sorted,'UniformOutput', false);
            bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           

            plot(time,bootstrap_summary(1,:),'Color', tt_colors{s,tt}); hold on;
            fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
                tt_colors{s,tt},'FaceAlpha', 0.3, 'Linestyle', 'none');
            plot(stimWindow, [24.5 24.5],'k','Linewidth',2)
            set(gca, 'box','off','TickDir','out')
            xlim(plotWindow);ylim([0 25]);
            ylabel('Firing rates (Hz)')
            xlabel('Time from stimulus onset (s)')
            title(titles{s})
        end   
    end
    sgtitle(recSites{i}) 
    % Save figure  
    zScoreFigPath = fullfile(mainDir,'Figures_grandAverage',[recSites{i} '_correct_trials']);
    print(zScoreFigPath,'-dpdf','-painters','-loose');
end

%% Peak responses
% Get peak activityof tHit and tCR during the stimulus window
mainDir = 'E:\Cross-Modal_Project\Revision\grand_average_and_peak_activity';

bin = 0.01;
timeWindow = [-1 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
stimWindow = [0 0.15];
stimOn_Bin = round(stimWindow(1)-timeWindow(1))*(1/bin);
stimOff_Bin = round((stimWindow(2)-timeWindow(1))*(1/bin)-1);
tts = {{1 1 1} {2 1 3}}; 
tt_names = {'tHit' 'tCR'};
tt_colors = {[0 0 1] [0 1 1]};

for i = 1: length(recSites)
    Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
    N = size(Xtrial,1);
    
    figure('Position',[0 0 800 400])
    for tt=1:length(tts)
        subplot(1,2,tt)
        Xtrial_tt = Xtrial(:,tts{tt}{1},tts{tt}{2},tts{tt}{3},stimOn_Bin:stimOff_Bin,:);
        X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials and different types of FA
        % Get peak response during the stimulus window
        X_tt = num2cell(X_tt,2); 
        peaks = cellfun(@(x) max(x,[],'all'), X_tt, 'UniformOutput', false);

        peak_hist = histogram(cell2mat(peaks),[0:10:160],'FaceColor',tt_colors{tt}, 'FaceAlpha',0.6);
        ylabel('Number of units')
        xlabel('Peak response (spike/sec)')
        xticks(0:10:160)
        xticklabels({0 [] [] [] 40 [] [] [] 80 [] [] [] 120 [] [] [] 160});
        set(gca, 'box','off','TickDir','out')
        title(tt_names{tt})
        if i == 4
            ylim([0 180])
            yticks(0:10:180)
            yticklabels({0 [] 20 [] 40 [] 60 [] 80 [] 100 [] 120 [] 140 [] 160 [] 180});
        else
            ylim([0 100])
            yticks(0:10:100)
            yticklabels({0 [] 20 [] 40 [] 60 [] 80 [] 100});
        end
    end   
    sgtitle(recSites{i}) 
    % Save figure  
    peakFigPath = fullfile(mainDir,'Figures_peaks',[recSites{i} '_peaks_v2']);
    print(peakFigPath,'-dpdf','-painters','-loose');
end

%% Grand average for all trial types (firing rates in Hz)
% 
% mainDir = 'E:\Cross-Modal_Project\Revision\grand_average_and_peak_activity';
% 
% bin = 0.01;
% timeWindow = [-1 2.5];
% time = timeWindow(1)+bin:bin:timeWindow(2);
% time_bw = [time, fliplr(time)];
% stimWindow = [0 0.15];
% nboot = 1000;
% plotWindow = [-1 2.5];
% 
% tts = {{1 1 1} {1 1 3} {2 1 1:2} {2 1 3};...
%     {2 2 2} {2 2 3} {1 2 1:2} {1 2 3};}; 
% 
% tt_names = {'Hit' 'Miss' 'FA' 'CR'};
% tt_colors= {[0 0 1] [0.15 0.15 0.15] [0.4 0.8 0.1] [1 0 0]};
% titles = {'Tactile stimulus' 'Visual stimulus'};
% 
% 
% for i = 1: length(recSites)
%     Xtrial = firingRates{strcmp(firingRates(:,1),recSites{i}),2};
%     N = size(Xtrial,1);
%     
%     figure('Position',[0 0 800 300])
%     for s=1:S
%         subplot(1,2,s)
%         for tt=1:length(tts)
%             Xtrial_tt = Xtrial(:,tts{s,tt}{1},tts{s,tt}{2},tts{s,tt}{3},:,:);
%             X_tt = squeeze(nanmean(Xtrial_tt,[2,3,4,6])); % average across trials and different types of FA
%             % bootstrapping for each time bin 
%             X_tt = num2cell(X_tt,1); 
%             bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), X_tt, 'UniformOutput', false);      
%             bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
%             bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)]', bootstat_sorted,'UniformOutput', false);
%             bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           
% 
%             plot(time,bootstrap_summary(1,:),'Color', tt_colors{tt}); hold on;
%             fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
%                 tt_colors{tt},'FaceAlpha', 0.3, 'Linestyle', 'none');
%             plot(stimWindow, [24.5 24.5],'k','Linewidth',2)
%             set(gca, 'box','off','TickDir','out')
%             xlim(plotWindow);ylim([0 25]);
%             ylabel('Firing rates (Hz)')
%             xlabel('Time from stimulus onset (s)')
%             title(titles{s})
%         end   
%     end
%     sgtitle(recSites{i}) 
%     % Save figure  
%     zScoreFigPath = fullfile(mainDir,'Figures_grandAverage',[recSites{i} '_all_trials_-1to2p5']);
%     print(zScoreFigPath,'-dpdf','-painters','-loose');
% end