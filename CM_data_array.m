% pre-stimulus lick and no lick trials
% 
% Construct data array
% Construct firing rate and trial number arrays from MSessionExplorer
% 1) Apply quality filters to select good units 
% 2) Remove outliers (mean firing rates > 50 Hz) 
% 3) Smooth the data (Gaussian kernel) for 10 ms bin
% 4) Remove pre-stimulus compulsive licking trials (lick within 1 sec before stimulus onsets)
% 5) Remove post-stimulus compulsive licking trials (lick within 0.1 sec after stimulus onsets)
% 6) Remove the sessions if the numbers of some trial types are less than 2
% except for TTV (respond to touch block, tactile stimulus, lick left) and
% VVT (respond to light block, visual stimulus, lick right)trials.
% 7) Remove sessions with poor performance (
    % 1. hit rate (visual and tactile stimuli) >= 35% 
    % 2. block performance (touch and light) >=  55% 
    % 3. overall performance >= 60%
    
% Outputs
% Data array for each session includes mouse id, date, recoding site, firing rate array, and trial number array.
% For each mouse, data array of sessions are concatenated and saved.

% MouseName - mouse id
% seshDate - session date
% recSite - recording site

% firingRates - firing rate is organized in a multidimensional array (neuron, stimulus, decision,
% context, time, trial) with size(N,B,S,D,T,K). Unit: Hz

% trialNum -  total trial number is organized in a multidimensional array (neuron, block, stimulus, decision) with size(N,B,S,D).

% earlyTrans - boolean array (1: early transition). For each block, the early transition is defined as 
% a window from the first trial of that block to the false alarm (included). It is organized in a multidimensional array 
% (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).

% lateTrans - boolean array (1: late transition). For each block, the late transition is defined as 
% a window from the first false alarm to the first hit(not included). It is organized in a multidimensional array 
% (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).

% beforeCue - boolean array (1: before a transition cue). It is organized in a multidimensional array 
% (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).

% afterCue -  boolean array (1: after a transition cue). It is organized in a multidimensional array 
% (neuron, block, stimulus, decision, trial) with size(N,B,S,D,K).

% behavPerformance - A table contains (1) the fractions of trial outcomes in respond-to-touch and respond-to-light blocks, 
% (2) the fractions of correct trials in respond-to-touch blocks, respond-to-light blocks, a session. 

% firstHit - the trial number of a first hit after block switch.  

% Note: A transition cue is a drop of water released from the rewarded port on the 9th trial of each block.
% It is not included in data array outputs. 

bin = 0.01;
timeWindow = [-1, 2.5];
T = (timeWindow(2) -timeWindow(1))/bin; % number of time points
B = 2; % number of block types: respond-to-touch and respond-to-light blocks
S = 2; % number of stimulus types: tactile and visual stimuli
D = 3; % number of decision types: right lick, left lick and no lick
L = 2; % number of pre-stimulus types: right and left lick
if bin == 0.01
    smoothWindow = 5; % smooth window of Gaussian filter 
end
preStimWindow = -1; % pre-stimulus onset no lick window
postStimWindow = 0.1; % post-stimulus onset no lick window % it was 0.15 before 090321
% preStimLickWindow = [-15 0]; % time window for detecting pre-stimulus licking
timeWindow_lick = [-0.5 0.5];
T_lick = (timeWindow_lick(2) -timeWindow_lick(1))/bin; % number of time points

load('E:\seFilePaths_SiProbe'); % File paths for MSessionExplorers
load('E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\filter'); % Unit quality filter

for m = 1:size(seFilePaths_all,1)
    tic
    MouseName = seFilePaths_all.MouseName{m};
    seFilePaths = seFilePaths_all.seFilePaths{m};
    % cluster ids for the MSessionExplorers of MouseName
    load(['E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\id_map\',MouseName,'_se_clusterIds']); 
    data_mouse = [];
    for i = 1: length(seFilePaths)
        load(seFilePaths{i})
        
        %%%%%% Load session and cluster information %%%%%%
        seshDate = se.userData.sessionInfo.seshDate;
        recSite = se.userData.sessionInfo.recSite;
        
        %%%%%% Get information about first hit %%%%%%
        summary_firstHit = CM.Process.GetFirstHitInfo(se);
        
        %%%%%% Pre-stimulus licking trials %%%%%%
        se_preStim_lick = se.Duplicate();
        rLickOnsetTimes = se_preStim_lick.GetColumn('behavTime', 'rLickOnset');
        lLickOnsetTimes = se_preStim_lick.GetColumn('behavTime', 'lLickOnset');
        isRlick = logical(cell2mat(cellfun(@(x) sum(x<0)>0, rLickOnsetTimes, 'UniformOutput',false)));
        isLlick = logical(cell2mat(cellfun(@(x) sum(x<0)>0, lLickOnsetTimes, 'UniformOutput',false)));
        isLick = isRlick | isLlick;
        se_preStim_lick.RemoveEpochs(~isLick);        
        % Align to the first lick during four second window preceeding stimulus
        rLickOnsetTimes_preStim = se_preStim_lick.GetColumn('behavTime', 'rLickOnset');
        lLickOnsetTimes_preStim = se_preStim_lick.GetColumn('behavTime', 'lLickOnset');
        isRlick_preStim = cell2mat(cellfun(@(x) sum(x<0)>0, rLickOnsetTimes_preStim, 'UniformOutput',false));
        isLlick_preStim = cell2mat(cellfun(@(x) sum(x<0)>0, lLickOnsetTimes_preStim, 'UniformOutput',false));
        [~,~,rfirstLick] = cellfun(@(x) find(x(x<0),1), rLickOnsetTimes_preStim, 'UniformOutput',false);
        [~,~,lfirstLick] = cellfun(@(x) find(x(x<0),1), lLickOnsetTimes_preStim, 'UniformOutput',false);
        ind_RLlick = find(sum([isRlick_preStim,isLlick_preStim],2)==2); % trials with both right and left licks
        if isempty(ind_RLlick)
            firstLick_preStim = cell2mat(cellfun(@(x,y) [x,y], rfirstLick, lfirstLick, 'UniformOutput',false));
        else
            for h=1:length(ind_RLlick) % trials with both right and left licks: use the first lick 
                if rfirstLick{ind_RLlick(h)} < lfirstLick{ind_RLlick(h)}
                    isLlick_preStim(ind_RLlick(h))= 0;
                    lfirstLick{ind_RLlick(h)} = [];
                else
                    isRlick_preStim(ind_RLlick(h))= 0;
                    rfirstLick{ind_RLlick(h)} = [];
                end    
            end
            firstLick_preStim = cell2mat(cellfun(@(x,y) [x,y], rfirstLick, lfirstLick, 'UniformOutput',false));
        end
        isLick_preStim = [isRlick_preStim isLlick_preStim];
        se_preStim_lick.AlignTime(firstLick_preStim)
        
        %%%%%% Remove trials with licking during the timewindow starting from -1 to 0.1 seconds, primarily due to compulsive licking %%%%%%        
%         rLickOnsetTimes = se.GetColumn('behavTime', 'rLickOnset');
%         lLickOnsetTimes = se.GetColumn('behavTime', 'lLickOnset');
        isRlick_comp = logical(cell2mat(cellfun(@(x) sum(x<postStimWindow & x>preStimWindow)>0, rLickOnsetTimes, 'UniformOutput',false)));
        isLlick_comp = logical(cell2mat(cellfun(@(x) sum(x<postStimWindow & x>preStimWindow)>0, lLickOnsetTimes, 'UniformOutput',false)));
        isLick_comp = isRlick_comp | isLlick_comp;        
        se.RemoveEpochs(isLick_comp);
        
        %%%%%% Get spike rates and remove units with poor quality and mean firing rate >50 (outlier) %%%%%% 
        cluster_ids = cluster_id_mouse{strcmp(cluster_id_mouse.seshDate, seshDate),3}{1};
        filter_session = filters(strcmp(filters.MouseName, MouseName) & strcmp(filters.seshDate, seshDate), :); 
        isPassed_session = [];
        for unit = 1: length(cluster_ids)
            isUnit = find(filter_session.cluster_id == cluster_ids(unit));
            isPassed_unit = filter_session.filter(isUnit);
            isPassed_session = [isPassed_session isPassed_unit];
        end  
        % Get spike rates
        spkRs = CM.Analysis.GetSpikeRate(se, bin, timeWindow(1), timeWindow(2)); 
        spkRs = spkRs(:,2:end); % remove time column
        spkRs_preStim_lick = CM.Analysis.GetSpikeRate(se_preStim_lick, bin, timeWindow_lick(1), timeWindow_lick(2)); 
        spkRs_preStim_lick = spkRs_preStim_lick(:,2:end); % remove time column
        % Apply quality filter and remove outliers
        isOutlier = mean(cell2mat(table2cell(spkRs)),1)*(1/bin)>50; % mean firing rates > 50 Hz
        spkRs = spkRs{:,logical(isPassed_session) & ~isOutlier};
        spkRs_preStim_lick = spkRs_preStim_lick{:,logical(isPassed_session) & ~isOutlier};
        
        %%%%%% Check number of trials for each trial type %%%%%%
        % Load trial type index
        trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
              'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
        ttInd = se.GetColumn('behavValue', trialTypes);
        tt_numb(1,:,:) = {1 2 3; 4 5 6}; % touch block- 2nd: stimulus types (tactile, visual), 3rd: decision types (lick right, lick left, no lick)
        tt_numb(2,:,:) = {7 8 9; 10 11 12}; % light block- 2nd: stimulus types (tactile, visual), 3rd: decision types (lick right, lick left, no lick)
        NumOfTrials = sum(ttInd,1);
        isMoreThan2Trials = sum(NumOfTrials([1,3,4,5,6,7,8,9,11,12])<2) == 0;
        
        ttInd_preStim_lick = se_preStim_lick.GetColumn('behavValue', trialTypes);
        
        %%%%%% Get Behavioral performance %%%%%%
        [isGoodPerf, behavPerformance] = CM.Process.GetBehavPerf(se);
        
        %%%%%% Get transition information %%%%%%
        % early and late transition 
        % before and after transition cue
        [isEarlyTrans, isLateTrans, isNotTrans, isBeforeCue, isAfterCue] = CM.Process.GetTransitionInfo(se);
        
        %%%%%% Get trial number information %%%%%%
        trialNumber = CM.Process.GetTrialNumber(se);
        
        %%%%%% Get lick onset %%%%%%
        firstLick = se.GetColumn('behavTime', 'firstLick');
        
        if ~isempty(spkRs) && isMoreThan2Trials && isGoodPerf
            N = size(spkRs,2); % number of neurons
            for n = 1:N
                for b = 1:B
                    for s=1:S
                        for d = 1:D
                            v = spkRs(ttInd(:,tt_numb{b,s,d}),n);
                            ind = find(ttInd(:,tt_numb{b,s,d}));
                            if ~isempty(v)
                                k = size(v,1); 
                                trialNum_total(n,b,s,d)= k; % total trial numbers 
                                for trial = 1:k
                                    % Convert to Hz
                                    if exist('smoothWindow') == 1 % smooth by Gaussian filter
                                        smoothFR = smoothdata(v{trial,1}*(1/bin),'gaussian',smoothWindow); 
                                        firingRates(n,b,s,d,1:T,trial) = [smoothFR];
                                    else
                                        FR = v{trial,1}*(1/bin);
                                        firingRates(n,b,s,d,1:T,trial) = [FR];
                                    end
                                    earlyTrans(n,b,s,d,trial) = isEarlyTrans(ind(trial));
                                    lateTrans(n,b,s,d,trial) = isLateTrans(ind(trial));
                                    notTrans(n,b,s,d,trial) = isNotTrans(ind(trial));
                                    beforeCue(n,b,s,d,trial) = isBeforeCue(ind(trial));
                                    afterCue(n,b,s,d,trial) = isAfterCue(ind(trial)); 
                                    lickOnset(n,b,s,d,trial) = firstLick(ind(trial)); 
                                    trialNumber_block(n,b,s,d,trial) = trialNumber(ind(trial)); 
                                end
                            end
                            
                            for l=1:L
                                v_preStim_lick = spkRs_preStim_lick(ttInd_preStim_lick(:,tt_numb{b,s,d}) & isLick_preStim(:,l),n);
                                if ~isempty(v_preStim_lick)
                                    k = size(v_preStim_lick,1); 
                                    trialNum_total_preStim_lick(n,b,s,d,l)= k; % total trial numbers 
                                    for trial = 1:k
                                        % Convert to Hz
                                        if exist('smoothWindow') == 1 % smooth by Gaussian filter
                                            smoothFR = smoothdata(v_preStim_lick{trial,1}*(1/bin),'gaussian',smoothWindow); 
                                            firingRates_preStim_lick(n,b,s,d,l,1:T_lick,trial) = [smoothFR];
                                        else
                                            FR = v_preStim_lick{trial,1}*(1/bin);
                                            firingRates_preStim_lick(n,b,s,d,l,1:T_lick,trial) = [FR];
                                        end                                 
                                    end
                                end
                            end
                        end
                    end
                end
            end
            new_data = [{MouseName} {seshDate} {recSite} {firingRates} {trialNum_total} {trialNumber_block}...
                {earlyTrans} {lateTrans} {notTrans} {beforeCue} {afterCue} {lickOnset} {behavPerformance} {summary_firstHit}...
                {firingRates_preStim_lick} {trialNum_total_preStim_lick}];
            data_mouse = [data_mouse; new_data];
        end
        clearvars -except seFilePaths_all seFilePaths data_mouse MouseName bin timeWindow T B S D time...
        filters cluster_id_mouse smoothWindow preStimWindow postStimWindow timeWindow_lick L T_lick % preStimLickWindow
    end
    if ~isempty(data_mouse)
        % Convert to a table
        VarNames = {'MouseName' 'seshDate' 'recSite' 'firingRate' 'Total_trialNum' 'trialNum_withinBlock' ...
            'earlyTrans' 'lateTrans' 'notTrans' 'beforeCue' 'afterCue' 'lickOnset' 'behavPerformance' 'firstHit'...
            'firingRate_preStim_lick' 'trialNum_preStim_lick'};
        data_mouse = cell2table(data_mouse,'VariableNames', VarNames); 

        % Save firing rate array for each mouse 
        if bin == 0.01
            dataPaths = fullfile('E:\Data_array\', ['AllTrials_', num2str(bin*1000), 'bin_Smooth_lickOnset'],...
                [MouseName,'_data_array']);
        else
            dataPaths = fullfile('E:\Data_array\', ['AllTrials_', num2str(bin*1000), 'bin_noSmooth_preStim_lick'],...
                [MouseName,'_data_array']);
        end
        save(dataPaths, 'data_mouse');
    end
    toc
end

%% Concatenate data array for all mice
dataDir = 'E:\Data_array\AllTrials_10bin_Smooth_lickOnset';
dataPaths = MBrowse.Files(dataDir);
data_all = table();
for i = 1: length(dataPaths)
    load(dataPaths{i})
    data_all = [data_all; data_mouse]; 
end
% Save zScore array for all mice 
dataSavePath = fullfile(dataDir, 'data_array');
save(dataSavePath, 'data_all','-v7.3');