% Raster plot of sensory-evoked activity
% Corrrect trial
% Chang et al., Figure 1f 

% Load data
load('E:\seFilePaths_SiProbe'); % File paths for MSessionExplorers
load('E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\filter'); % Unit quality filter
load('E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity\Unit_AUC_tCorrect') % tHit and tCR selectivity

% Setting
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};

bin = 0.01;
timeWindow = [-1, 2.5];
time = timeWindow(1)+bin:bin:timeWindow(2);
time_bw = [time, fliplr(time)];
smoothWindow = 5; % smooth window of Gaussian filter 
preStimWindow = -1; % pre-stimulus onset no lick window
postStimWindow = 0.1; % post-stimulus onset no lick window 
tt_colors = {[0 0 1] [0 1 1]; [1 0 0] [1 0 1]};
block_bar_colors = {[0.05 0.5 0.1] [0.9 0.8 0.2];[0.9 0.8 0.2] [0.05 0.5 0.1]};
plotWindow = [-0.25 1];
title_names = {'Tactile stimulus' 'Visual stimulus'};
stimWindow = [0 0.15];
nboot=100;
save_path = 'E:\CM_NeuralActivity_Analysis\RasterPlot\Response_correctTrials';

for i=1: length(recSites)
    AUC_recSite = AUC_tCorrect(strcmp(AUC_tCorrect(:,3), recSites{i}),:);
    recSite = recSite_names{i};
    activity_recSite = [];
    Unit_number = 0;
    for session = 1:length(AUC_recSite)
        MouseName = AUC_recSite{session,1};
        seshDate = AUC_recSite{session, 2};
        AUC = AUC_recSite{session, 4};
        
        % Find corresponding MSessionExplorers 
        seFilePaths = seFilePaths_all{strcmp(seFilePaths_all.MouseName,MouseName),3}{1,1};
        isSession = cell2mat(cellfun(@(x) strcmp(x(end-12:end-7),seshDate), seFilePaths, 'UniformOutput',false));
        load(seFilePaths{isSession})
        % cluster ids for the MSessionExplorers of MouseName
        load(['E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\id_map\',MouseName,'_se_clusterIds']); 
        
        %%%%%% Remove trials with licking during the timewindow starting from -1 to 0.1 seconds, primarily due to compulsive licking %%%%%%        
        rLickOnsetTimes = se.GetColumn('behavTime', 'rLickOnset');
        lLickOnsetTimes = se.GetColumn('behavTime', 'lLickOnset');
        isRlick_comp = logical(cell2mat(cellfun(@(x) sum(x<postStimWindow & x>preStimWindow)>0, rLickOnsetTimes, 'UniformOutput',false)));
        isLlick_comp = logical(cell2mat(cellfun(@(x) sum(x<postStimWindow & x>preStimWindow)>0, lLickOnsetTimes, 'UniformOutput',false)));
        isLick_comp = isRlick_comp | isLlick_comp;
        se.RemoveEpochs(isLick_comp);
        
        %%%%%% Get spike rates and remove units with poor quality and mean firing rate >50 (outlier) %%%%%% 
        cluster_ids = cluster_id_mouse{strcmp(cluster_id_mouse.seshDate, seshDate),3}{1};
        filter_session = filters(strcmp(filters.MouseName, MouseName) & strcmp(filters.seshDate, seshDate), :); 
        isPassed_session = [];
        for n = 1: length(cluster_ids)
            isUnit = find(filter_session.cluster_id == cluster_ids(n));
            isPassed_unit = filter_session.filter(isUnit);
            isPassed_session = [isPassed_session isPassed_unit];
        end  
        % Get spike and spike rate data
        spk = se.GetTable('spikeTime');
        spkRs = CM.Analysis.GetSpikeRate(se, bin, timeWindow(1), timeWindow(2)); 
        spkRs = spkRs(:,2:end); % remove time column
        % Apply quality filter and remove outliers
        isOutlier = mean(cell2mat(table2cell(spkRs)),1)*(1/bin)>50; % mean firing rates > 50 Hz
        spk = spk{:,logical(isPassed_session) & ~isOutlier};  
        spkRs = spkRs{:,logical(isPassed_session) & ~isOutlier};
        N_session = size(spk,2); % number of neurons
        
        % Load trial type index
        trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                     'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
        ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
                'VariableNames', trialTypes);
        new_ttInd = {ttInd.TTT ttInd.VTN;
            ttInd.VVV ttInd.TVN};  % row: tactile and visual stimulus; column: hit, cr
        
        % reaction time 
        lickOn = se.GetColumn('behavTime', 'firstLick');

        % Check if the number of units of MSessionExplorer == the number of units of AUC data
        if N_session == length(AUC)
            for n = 1:N_session 
                s = spk(:,n); % Select unit
                sr = spkRs(:,n);
                figure('Position', [0,0, 480, 600]);clf % set figure size
                % Spike raster plots and lick onsets
                for stim = 1:size(new_ttInd,1)
                    switch stim
                        case 1
                            subplot('Position', [0.2 0.35 0.25 0.5])
                            Hit_color = 'b';
                        case 2
                            subplot('Position', [0.6 0.35 0.25 0.5])
                            Hit_color = 'r';
                    end
                    base=0;
                    for tt = 1:size(new_ttInd,2)
                        s_isTrial = s(new_ttInd{stim,tt});
                        if tt==1
                            lickOn_isTrial = lickOn(new_ttInd{stim,tt});         
                            [lickOn_isTrial,sort_idx] = sortrows(lickOn_isTrial);
                            s_isTrial = s_isTrial(sort_idx);
                        end
                        ind = 1+base:size(s_isTrial,1)+base;
                        MPlot.PlotRaster(s_isTrial, ind, 0.5, 'Color', 'k');
                        if tt==1
                            scatter(lickOn_isTrial,ind, 8, 'filled', Hit_color);
                        end
                        base = base + size(s_isTrial,1);
                        trialNum_summary(stim,tt) = length(s_isTrial);
                    end
                end
                tStim_trialNum = sum(trialNum_summary(1,:));
                vStim_trialNum = sum(trialNum_summary(2,:));
                % block length
                [M,I] = max([tStim_trialNum vStim_trialNum]);
                if I == 1
                    tStim_barRatio = tStim_trialNum/tStim_trialNum;
                    vStim_barRatio = vStim_trialNum/tStim_trialNum;
                else
                    tStim_barRatio = tStim_trialNum/vStim_trialNum;
                    vStim_barRatio = vStim_trialNum/vStim_trialNum;
                end

                % figure setting & block type bar
                for stim = 1:size(new_ttInd,1)
                    % setting of raster plots
                    switch stim
                        case 1
                            subplot('Position', [0.2 0.35 0.25 0.5])
                        case 2
                            subplot('Position', [0.6 0.35 0.25 0.5])
                    end
                    xlim(plotWindow);
                    ylim([0.5 M+3]);
                    xticks([-0.2:0.2:1])
                    ylabel('Trials');
                    title(title_names{stim}, 'FontSize',8);
                    set(gca, 'box','off','TickDir','out','xticklabel',{[]})
                    plot(stimWindow, [M+1 M+1],'k','Linewidth',2)

                    % block type bar
                    switch stim
                        case 1
                            subplot('Position', [0.453 0.35 0.05 0.5*tStim_barRatio])
                        case 2
                            subplot('Position', [0.853 0.35 0.05 0.5*vStim_barRatio])
                    end
                    x1 = repmat(0,[sum(trialNum_summary(stim,1)),1]); % hit
                    y1 = (1:sum(trialNum_summary(stim,1))).';
                    x2 = repmat(0,[sum(trialNum_summary(stim,2)),1]); % cr
                    y2 = (sum(trialNum_summary(stim,1))+1:sum(trialNum_summary(stim,:))).';
                    scatter(x1,y1,12,block_bar_colors{stim,1},'s','filled');
                    hold on
                    scatter(x2,y2,12,block_bar_colors{stim,2},'s','filled');
                    ylim([0.5 sum(trialNum_summary(stim,:))+0.5]);
                    set(gca,'visible','off')
                end

                % Average trace
                for stim = 1:size(new_ttInd,1)
                    switch stim
                        case 1
                            subplot('Position', [0.2 0.1 0.25 0.2])
                        case 2
                            subplot('Position', [0.6 0.1 0.25 0.2])
                    end

                    for tt = 1:size(new_ttInd,2)
                        sr_isTrial = sr(new_ttInd{stim,tt});
    %                     sr_isTrial_hz = cellfun(@(x) (x*(1/bin)).', sr_isTrial,'UniformOutput', false);
    %                     sr_isTrial_reshpaed = num2cell(cell2mat(sr_isTrial_hz),1);
                        sr_isTrial_smoothed = cellfun(@(x) smoothdata(x*(1/bin),'gaussian',smoothWindow).', sr_isTrial,'UniformOutput', false);
                        sr_isTrial_reshpaed = num2cell(cell2mat(sr_isTrial_smoothed),1);
                        % bootstrapping for each time bin
                        bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), sr_isTrial_reshpaed, 'UniformOutput', false);      
                        bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
                        bootstrap_summary = cellfun(@(x) [mean(x), x(round(0.025*nboot)), x(round(0.975*nboot))]', bootstat_sorted,'UniformOutput', false);
                        bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI  
                        plot(time,bootstrap_summary(1,:),'Color', tt_colors{stim,tt}); hold on;
                        fill(time_bw, [bootstrap_summary(2,:), fliplr(bootstrap_summary(3,:))],...
                            tt_colors{stim,tt},'FaceAlpha', 0.3, 'Linestyle', 'none');
                        activity_unit(stim,tt) = {bootstrap_summary};
                    end
                end
                for stim = 1:size(new_ttInd,1)
                    switch stim
                        case 1
                            subplot('Position', [0.2 0.1 0.25 0.2])
                        case 2
                            subplot('Position', [0.6 0.1 0.25 0.2])
                    end
                    set(gca, 'box','off','TickDir','out')
                    xlim(plotWindow);
                    ylim_max = max(cell2mat(activity_unit),[],'all');
                    ylim([0 ceil((ylim_max+0.1*ylim_max))]);
                    xticks([-0.2:0.2:1])
                    ylabel('Spikes/s')
                    xlabel('Time from stimulus onset (s)')
                end    

                sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(n+Unit_number)],'FontSize',10);
                activity_recSite = [activity_recSite;{activity_unit}];
                % Save fig
                FigPath = fullfile(save_path, recSite,[MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(n+Unit_number), '.pdf']);
                print(FigPath,'-dpdf','-painters','-loose');
            end
            Unit_number = Unit_number +N_session;
            close all
        end
    end
    activity{i,1} = recSite;
    activity{i,2} = activity_recSite;
end

