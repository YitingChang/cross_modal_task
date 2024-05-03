% Raster plot of baseline activity 
% Chang et al., Fig 2a

% Load data
load('E:\CM_NeuralActivity_Analysis\ROC\block_selectivity\Unit_AUC_BP_100msBin_1000nBoot') % block selectivity
load('E:\seFilePaths_SiProbe'); % File paths for MSessionExplorers
load('E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\filter'); % Unit quality filter

% Setting
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
bin = 0.1;
timeWindow = [-1, 2.5];
B = 2; % number of blocks: respond-to-touch and respond-to-light blocks
S = 2; % number of stimuli: tactile and visual stimuli
D = 3; % number of decisions: right lick, left lick and no lick
preStimWindow = -1; % pre-stimulus onset no lick window
postStimWindow = 0.1; % post-stimulus onset no lick window 
initialBin = (0-timeWindow(1))*(1/bin);
nBoot = 1000;
alpha = 0.05;
save_path = 'E:\CM_NeuralActivity_Analysis\RasterPlot\Baseline_allTrials';
for i=1: length(recSites)
    AUC_recSite = AUC_rule(strcmp(AUC_rule(:,3), recSites{i}),:);
    N = sum(cell2mat(cellfun(@(x) size(x,2),AUC_recSite(:,4),'UniformOutput',false))); 
    alpha_Bonferroni = alpha/N; % Bonferroni correction
    recSite = recSite_names{i};
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
        N_session = size(spk,2); % number of neurons
        
        % Create index for block numbers
        istBlock = ismember(se.GetColumn('behavValue', 'blockType'),'Whisker');
        isvBlock = ismember(se.GetColumn('behavValue', 'blockType'),'Visual');
        isChange = ismember(se.GetColumn('behavValue', 'blockType'), 'Change') ;

        % Check if the number of units of MSessionExplorer == the number of units of AUC data
        if N_session == length(AUC)
            for n = 1:N_session 
                %%% AUC %%%
                AUC_unit = AUC{n};
                AUC_mean = round(mean(AUC_unit(initialBin,:)),2);
                AUC_upCI = AUC_unit(initialBin,ceil((1-alpha_Bonferroni/2)*nBoot));
                AUC_lowCI = AUC_unit(initialBin,ceil((alpha_Bonferroni/2)*nBoot));
                isSignificant = (0.5-AUC_upCI)*(0.5-AUC_lowCI)>0;

                %%% Plotting %%%
                % Raster plot
                s = spk(:,n); % Select unit
                figure('Position', [0,0, 200, 800]); % set figure size
                subplot('Position', [0.25 0.1 0.55 0.8])
                ind = 1:size(s,1);
                MPlot.PlotRaster(s(ind), ind, 0.5, 'Color', 'k'); hold on
                xlim([-1 0]);
                ylim([ind(1)-0.5 length(s)+0.5]);
                xlabel('Time from stimulus onset');
                ylabel('Trial number within a session');
                set(gca, 'box','off','TickDir','out')
                
                % Block type bar
                subplot('Position', [0.845 0.1 0.05 0.8])
                x1 = repmat(0,[sum(istBlock),1]);
                y1 = find(istBlock == 1);
                x2 = repmat(0,[sum(isvBlock),1]);
                y2 = find(isvBlock == 1);
                x3 = repmat(-0.5,[sum(isChange),1]);
                y3 = find(isChange ==1);
                scatter(x1,y1,12,'b','s','filled');hold on
                scatter(x2,y2,12,'r','s','filled');
                scatter(x3,y3,24,'k','<','filled');
                ylim([ind(1)-0.5 length(s)+0.5]);
                set(gca,'visible','off')
                sgtitle({[MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(Unit_number+n)],...
                    ['AUC=', num2str(AUC_mean), '\_sig\_',num2str(isSignificant)]},'FontSize',10);

                % Save fig
                rasterFigPath = fullfile(save_path,recSite,[MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(Unit_number+n),...
                    '_sig_',num2str(isSignificant),'.pdf']);
                print(rasterFigPath,'-dpdf','-painters','-loose');
            end
        else
            disp(['Error: ' MouseName ' ' seshDate ' ' recSite])
        end
        Unit_number = Unit_number + N_session;
    end
    close all
end