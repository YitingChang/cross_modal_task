%% tHit vs. tCR
% Touch-evoked responses
% Chang et al., Fig 1h

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% Import setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
nBoot = 1000;
rng(15)
pos_tt = {[1 1 1]}; % tHit
neg_tt = {[2 1 3]}; % tCR
AUC_tCorrect = [];
for i = 1: length(recSites)
    u = data_all(strcmp(data_all.recSite, recSites{i}), :);
    AUC_recSite = u{:,1:3};
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1};
        N = size(firingRates_session,1); % number of neurons
        trialNum_session = u{session,5}{1,1};
        parfor n = 1:N % parallel for loop 
            firingRates_pos = [];
            firingRates_neg = [];
            for tt=1:length(pos_tt)
                firingRates_pos_tt = squeeze(firingRates_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3),:,...
                    1:trialNum_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3))));
                firingRates_neg_tt = squeeze(firingRates_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3),:,...
                    1:trialNum_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3))));
                firingRates_pos = [firingRates_pos firingRates_pos_tt];
                firingRates_neg = [firingRates_neg firingRates_neg_tt];
            end  
            
            posTrials = num2cell(firingRates_pos,2);
            negTrials = num2cell(firingRates_neg,2);
            posTrialNum = size(firingRates_pos,2);
            negTrialNum = size(firingRates_neg,2);
            binNum = size(posTrials,1);
            % bootstrapped posTrials and negTrials
            [~, posTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), posTrials, 'UniformOutput', false); % Get bootstrapped indexes 
            posTrials_bootstrapped = cellfun(@(x,y) x(y), posTrials, posTrials_bootstrapped_ind, 'UniformOutput', false); % Get bootstrapped values
            posTrials_bootstrapped_matrix = mat2cell(cell2mat(posTrials_bootstrapped),[repmat(posTrialNum,[binNum 1])], [ones([nBoot,1])]); % reshape posTrials_bootstrapped for AUC calculation
            [~, negTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), negTrials, 'UniformOutput', false);  
            negTrials_bootstrapped = cellfun(@(x,y) x(y), negTrials, negTrials_bootstrapped_ind, 'UniformOutput', false);
            negTrials_bootstrapped_matrix = mat2cell(cell2mat(negTrials_bootstrapped),[repmat(negTrialNum,[binNum 1])], [ones([nBoot,1])]); 

            Labels = [ones(posTrialNum,1); zeros(negTrialNum,1)];
            [~,~,~,AUC] = cellfun(@(x,y) perfcurve(Labels,[x(:); y(:)],1),...
                               posTrials_bootstrapped_matrix, negTrials_bootstrapped_matrix, 'uni', 0);
                           
            AUC_sorted = sort(cell2mat(AUC),2); 
            AUC_session{n} = AUC_sorted; % AUC: timeBin x nBoot
        end
        AUC_recSite{session,4} = AUC_session;
        clear AUC_session
    end
    AUC_tCorrect = [AUC_tCorrect; AUC_recSite];
    clear AUC_recSite
end
aucPaths = fullfile(mainDir,['Unit_AUC_tCorrect_' num2str(nBoot) 'nBoot']);
save(aucPaths, 'AUC_tCorrect');

%% histogram
% distribution of mean AUC of touch-evoked responses
% label significant units
% normalization: probability
bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow = [0 0.15];
binStart =round((analysisWindow(1)-timeWindow(1))/bin)+1;
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);
save_path= ('E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity');
figure('Position',[0 0 200 800])
nBoot =1000;
load('E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity\Unit_AUC_tCorrect_1000nBoot')
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
alpha = 0.05;

for i=1: length(recSites)
    AUC_recSite = AUC_tCorrect(strcmp(AUC_tCorrect(:,3), recSites{i}), 4);
    AUC_sessions = [];
    % concatenation across sessions
    for session=1:size(AUC_recSite,1)
        AUC_session = AUC_recSite{session};
        AUC_sessions = [AUC_sessions AUC_session];
    end
    N = size(AUC_sessions,2);
    alpha_Bonferroni = alpha/N; % Bonferroni correction
    AUC_mean = cellfun(@(x) mean(x(binStart:binEnd,:),'all'), AUC_sessions,'UniformOutput', false);
    AUC_upCI = cellfun(@(x) mean(x(binStart:binEnd,ceil((1-alpha_Bonferroni/2)*nBoot))), AUC_sessions,'UniformOutput', false);
    AUC_lowCI = cellfun(@(x) mean(x(binStart:binEnd,ceil((alpha_Bonferroni/2)*nBoot))), AUC_sessions,'UniformOutput', false);
    isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
    sig_percent = round(sum(isSignificant)/length(isSignificant)*100,1); 
    
    subplot(4,1,i)
    h_all = histogram(cell2mat(AUC_mean),[0:0.05:1], 'FaceColor','white');hold on;
    h_sig = histogram(cell2mat(AUC_mean(isSignificant)),[0:0.05:1], 'FaceColor','black', 'FaceAlpha',0.6);hold off; 
    h_all_normalization = h_all.Values/length(AUC_mean)*100;
    h_sig_normalization = h_sig.Values/length(AUC_mean)*100;
    h_all_prb = histogram('BinEdges',[0:0.05:1], 'BinCounts',h_all_normalization,'FaceColor','white');hold on;
    h_sig_prb = histogram('BinEdges',[0:0.05:1], 'BinCounts',h_sig_normalization,'FaceColor','black', 'FaceAlpha',0.6);hold off; 
    ylim([0 40])
    ylabel('Percentage of units')
    yticks(0:5:40)
    yticklabels({0 [] 10 [] 20 [] 30 [] 40});
    xlim([0 1])
    xlabel('AUC')
    xticks(0:0.1:1)
    xticklabels({0 [] [] [] [] 0.5 [] [] [] [] 1});
    set(gca, 'box','off','TickDir','out')
    text(0.75,25,[num2str(sig_percent), '%'])
    title(recSites{i})
end

% Save fig
histFigPath = fullfile(save_path, ['AUC_histogram_tCorrect_Bonferroni_' num2str(nBoot) 'nBoot_p0top15.pdf']);
print(histFigPath,'-dpdf','-painters','-loose');

%% Heatmap for tHit-tCR selectivity across time
bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow = [0 0.25];
binStart =(analysisWindow(1)-timeWindow(1))/bin+1;
binEnd = (analysisWindow(2)-timeWindow(1))/bin;
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
% xTickLables = [repmat(' ',99,1);'0'; repmat(' ',99,1); '1'; repmat(' ',99,1); '2'; repmat(' ',50,1)]; 
xTickLables = ['0'; repmat(' ',9,1); '1'; repmat(' ',9,1); '2'; repmat(' ',4,1)]; 

alpha = 0.05;
nBoot =1000;

for b= binStart:binEnd
    for i=1: length(recSites)
        AUC_recSite = AUC_tCorrect(strcmp(AUC_tCorrect(:,3), recSites{i}), 4);
        AUC_sessions = [];
        % concatenation across sessions
        for session=1:size(AUC_recSite,1)
            AUC_session = AUC_recSite{session};
            AUC_sessions = [AUC_sessions AUC_session];
        end
        N = size(AUC_sessions,2);
        alpha_Bonferroni = alpha/N; % Bonferroni correction
        AUC_mean = cellfun(@(x) mean(x(b,:)), AUC_sessions,'UniformOutput', false);
        AUC_upCI = cellfun(@(x) x(b,ceil((1-alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        AUC_lowCI = cellfun(@(x) x(b,ceil((alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
        sig_percent = sum(isSignificant)/length(isSignificant);
        probability_heatmap(i,b-binStart+1) = sig_percent;      
    end
end

figure('Position', [0 0 600 300])
h= heatmap(probability_heatmap, 'Colormap',bone, 'ColorLimits',[0 0.5],...%'GridVisible','off',...
                'ColorbarVisible', 'on','Title','tHit vs tCR selectivity');
h.XDisplayLabels = xTickLables; 
h.YDisplayLabels = recSite_names;
h.XLabel = 'Time from stimulus onset (s)';

% Save fig
save_path= ('E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity'); 
FigPath = fullfile(save_path, ['tHit_tCR_selectivity_Bonferroni_p0top25' num2str(nBoot) 'nBoot.pdf']);
print(FigPath,'-dpdf','-painters','-loose');
%% Correct trials
% tBlock vs. vBlock
% Baseline activity of correct trials in respond-to-touch block vs. respond-to-light blocks
% Chang et al., Fig 3c

% tStim vs. vStim
% Negative controls: baseline activity of correct tactile trials vs. visual trials 
% Positive controls: sensory-evoked activity of correct tactile trials vs. visual trials 
% Chang et al., Fig. S3

% load data array
load('E:\Data_array\AllTrials_100bin_noSmooth\data_array.mat')
% Import setting 
% mainDir = 'E:\CM_NeuralActivity_Analysis\ROC\block_selectivity';
mainDir = 'E:\CM_NeuralActivity_Analysis\ROC\stimulus_selectivity';

bin = 0.1;
% component_names = {'stimulus' 'decision' 'context'};
nBoot = 1000;
rng(15)
timeWindow = [-1 2.5];
analysisWindow = [-1 0.3];
binStart =round((analysisWindow(1)-timeWindow(1))/bin+1);
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);

recSites = {'left S1', 'left S2', 'left wM2','left ALM'};

% block probability
pos_tt = {[1 1 1] [1 2 3]}; % tHit, vCR
neg_tt = {[2 1 3] [2 2 2]}; % tCR, vHit

% stimulus probability
% pos_tt = {[1 1 1] [2 1 3]}; % tHit, tCR
% neg_tt = {[2 2 2] [1 2 3]}; % vHit, vCR

% decision probability (merged left and right lick)
% pos_tt = {[1 1 1] [2 2 2]}; % tHit, vHit
% neg_tt = {[1 3 2] [2 3 1]}; % tCR, vCR

AUC_stimulus = [];
for i = 1: length(recSites)
    u = data_all(strcmp(data_all.recSite, recSites{i}), :);
    AUC_recSite = u{:,1:3};
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1}(:,:,:,:,binStart:binEnd,:);
        N = size(firingRates_session,1); % number of neurons
        trialNum_session = u{session,5}{1,1};
        parfor n = 1:N % parallel for loop 
            firingRates_pos = [];
            firingRates_neg = [];
            for tt=1:length(pos_tt)
                firingRates_pos_tt = squeeze(firingRates_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3),:,...
                    1:trialNum_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3))));
                firingRates_neg_tt = squeeze(firingRates_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3),:,...
                    1:trialNum_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3))));
                firingRates_pos = [firingRates_pos firingRates_pos_tt];
                firingRates_neg = [firingRates_neg firingRates_neg_tt];
            end  
            
            posTrials = num2cell(firingRates_pos,2);
            negTrials = num2cell(firingRates_neg,2);
            posTrialNum = size(firingRates_pos,2);
            negTrialNum = size(firingRates_neg,2);
            binNum = size(posTrials,1);
            % bootstrapped posTrials and negTrials
            [~, posTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), posTrials, 'UniformOutput', false); % Get bootstrapped indexes 
            posTrials_bootstrapped = cellfun(@(x,y) x(y), posTrials, posTrials_bootstrapped_ind, 'UniformOutput', false); % Get bootstrapped values
            posTrials_bootstrapped_matrix = mat2cell(cell2mat(posTrials_bootstrapped),[repmat(posTrialNum,[binNum 1])], [ones([nBoot,1])]); % reshape posTrials_bootstrapped for AUC calculation
            [~, negTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), negTrials, 'UniformOutput', false);  
            negTrials_bootstrapped = cellfun(@(x,y) x(y), negTrials, negTrials_bootstrapped_ind, 'UniformOutput', false);
            negTrials_bootstrapped_matrix = mat2cell(cell2mat(negTrials_bootstrapped),[repmat(negTrialNum,[binNum 1])], [ones([nBoot,1])]); 

            Labels = [ones(posTrialNum,1); zeros(negTrialNum,1)];
            [~,~,~,AUC] = cellfun(@(x,y) perfcurve(Labels,[x(:); y(:)],1),...
                               posTrials_bootstrapped_matrix, negTrials_bootstrapped_matrix, 'uni', 0);
                           
            AUC_sorted = sort(cell2mat(AUC),2);
            AUC_session{n} = AUC_sorted; % AUC: timeBin x nBoot
        end
        AUC_recSite{session,4} = AUC_session;
        clear AUC_session
    end
    AUC_stimulus = [AUC_stimulus; AUC_recSite];
    clear AUC_recSite
end
aucPaths = fullfile(mainDir,['Unit_AUC_BP_100msBin_' num2str(nBoot) 'nBoot_-1top3']);
save(aucPaths, 'AUC_block');
 

%% histogram
% distribution of AUC 
% label significant units
% normalization: probability
bin = 0.1; % sec
analysisWindow = [-1 2.5];
binNum = (analysisWindow(2) - analysisWindow(1))/bin;
nBoot =1000;
% load('E:\CM_NeuralActivity_Analysis\ROC\block_selectivity\Unit_AUC_BP_100msBin_10000nBoot_-p3top3')
load('E:\CM_NeuralActivity_Analysis\ROC\block_selectivity\Unit_AUC_BP_100msBin_1000nBoot_-1to2p5')
% load('E:\CM_NeuralActivity_Analysis\ROC\stimulus_selectivity\Unit_AUC_SP_100msBin_1000nBoot_-p3top3')

recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
alpha = 0.05;

for b=1:binNum
    figure('Position',[0 0 800 200])
    for i=1: length(recSites)
        AUC_recSite = AUC_rule(strcmp(AUC_rule(:,3), recSites{i}), 4);
        AUC_sessions = [];
        session_id_all = [];
        % concatenation across sessions
        for session=1:size(AUC_recSite,1)
            AUC_session = AUC_recSite{session};
            AUC_sessions = [AUC_sessions AUC_session];
            n = size(AUC_session,2);
            session_id = repmat(session, [n,1]);
            session_id_all = [session_id_all; session_id];
        end
        N = size(AUC_sessions,2);
        alpha_Bonferroni = alpha/N; % Bonferroni correction
        AUC_mean = cellfun(@(x) mean(x(b,:)), AUC_sessions,'UniformOutput', false);
        AUC_upCI = cellfun(@(x) x(b,ceil((1-alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        AUC_lowCI = cellfun(@(x) x(b,ceil((alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
%         test = session_id_all(isSignificant);
        sig_percent = round(sum(isSignificant)/length(isSignificant)*100,1);
        subplot(1,4,i)
        h_all = histogram(cell2mat(AUC_mean),[0:0.05:1], 'FaceColor','white');hold on;
        h_sig = histogram(cell2mat(AUC_mean(isSignificant)),[0:0.05:1], 'FaceColor','black', 'FaceAlpha',0.6);hold off; 
        h_all_normalization = h_all.Values/length(AUC_mean);
        h_sig_normalization = h_sig.Values/length(AUC_mean);
        h_all_prb = histogram('BinEdges',[0:0.05:1], 'BinCounts',h_all_normalization,'FaceColor','white');hold on;
        h_sig_prb = histogram('BinEdges',[0:0.05:1], 'BinCounts',h_sig_normalization,'FaceColor','black', 'FaceAlpha',0.6);hold off; 
        max_count = max([h_all_prb.Values h_sig_prb.Values]);
        min_count = min([h_all_prb.Values h_sig_prb.Values]);
        ylim_max = round(max_count+ 0.1*(max_count - min_count),2);
    %     ylim([0 ylim_max]);
        ylim([0 0.5]);
        xlabel('AUC')
        xticks([0:0.1:1])
        xticklabels({0 [] [] [] [] 0.5 [] [] [] [] 1})
    %     yticks([0:0.05:ylim_max])
        yticks([0:0.1:0.6])
        set(gca, 'box','off','TickDir','out')
        text(0.75,0.3,[num2str(sig_percent), '%'])
        if i==1
            ylabel('Proportion of units')
        end
        title(recSites{i})
        
        block_preference{i,1} = sum(h_sig_normalization(11:end));% tBlock 
        block_preference{i,2} = sum(h_sig_normalization(1:10));% vBlock
    end
    startTime = round(analysisWindow(1) + (b-1)*bin,1); 
    endTime =  round(analysisWindow(1) + b*bin,1);
    sgtitle([num2str(startTime) ' ~ ' num2str(endTime) ' s']);

    % Save fig
    save_path= ('E:\CM_NeuralActivity_Analysis\ROC\block_selectivity\BP_histogram_100msBin_1000nBoot_Bonferroni_-1to2p5'); 
%     histFigPath = fullfile(save_path,['BP_histogram_' num2str(startTime) 'to' num2str(endTime) 's.pdf']);
%     save_path= ('E:\CM_NeuralActivity_Analysis\ROC\stimulus_selectivity\SP_histogram_100msBin_1000nBoot_Bonferroni'); 
    histFigPath = fullfile(save_path,['BP_histogram_' num2str(startTime) 'to' num2str(endTime) 's.pdf']);
    print(histFigPath,'-dpdf','-painters','-loose');
end


%% stimulus probability: CR vs baseline or nolick vs baseline

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% Import setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\ROC\10ms_bin\stimulus_probability';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
nBoot = 1000;
rng(15)

bin = 0.01; % sec
timeWindow = [-1 2.5];
posWindow = [0 0.25];
negWindow = [-0.25 0];
posStart =(posWindow(1)-timeWindow(1))/bin+1;
posEnd = (posWindow(2)-timeWindow(1))/bin;
negStart =(negWindow(1)-timeWindow(1))/bin+1;
negEnd = (negWindow(2)-timeWindow(1))/bin;

% SP for tactile trials
% pos_tt = {[2 1 3] [1 1 3]}; % no lick
% neg_tt = {[2 1 3] [1 1 3]}; % no lick baseline
% SP for visual trials
pos_tt = {[2 2 3] [1 2 3]}; % no lick
neg_tt = {[2 2 3] [1 2 3]}; % no lick baseline

AUC_vSP = [];
for i = 1: length(recSites)
    u = data_all(strcmp(data_all.recSite, recSites{i}), :);
    AUC_recSite = u{:,1:3};
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1};
        N = size(firingRates_session,1); % number of neurons
        trialNum_session = u{session,5}{1,1};
        parfor n = 1:N % parallel for loop 
            firingRates_pos = [];
            firingRates_neg = [];
            for tt=1:length(pos_tt)
                firingRates_pos_tt = squeeze(firingRates_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3),posStart:posEnd,...
                    1:trialNum_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3))));
                firingRates_neg_tt = squeeze(firingRates_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3),negStart:negEnd,...
                    1:trialNum_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3))));
                firingRates_pos = [firingRates_pos firingRates_pos_tt];
                firingRates_neg = [firingRates_neg firingRates_neg_tt];
            end  
            
            posTrials = num2cell(firingRates_pos,2);
            negTrials = num2cell(firingRates_neg,2);
            posTrialNum = size(firingRates_pos,2);
            negTrialNum = size(firingRates_neg,2);
            binNum = size(posTrials,1);
            % bootstrapped posTrials and negTrials
            [~, posTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), posTrials, 'UniformOutput', false); % Get bootstrapped indexes 
            posTrials_bootstrapped = cellfun(@(x,y) x(y), posTrials, posTrials_bootstrapped_ind, 'UniformOutput', false); % Get bootstrapped values
            posTrials_bootstrapped_matrix = mat2cell(cell2mat(posTrials_bootstrapped),[repmat(posTrialNum,[binNum 1])], [ones([nBoot,1])]); % reshape posTrials_bootstrapped for AUC calculation
            [~, negTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), negTrials, 'UniformOutput', false);  
            negTrials_bootstrapped = cellfun(@(x,y) x(y), negTrials, negTrials_bootstrapped_ind, 'UniformOutput', false);
            negTrials_bootstrapped_matrix = mat2cell(cell2mat(negTrials_bootstrapped),[repmat(negTrialNum,[binNum 1])], [ones([nBoot,1])]); 

            Labels = [ones(posTrialNum,1); zeros(negTrialNum,1)];
            [~,~,~,AUC] = cellfun(@(x,y) perfcurve(Labels,[x(:); y(:)],1),...
                               posTrials_bootstrapped_matrix, negTrials_bootstrapped_matrix, 'uni', 0);
                           
            AUC_sorted = sort(cell2mat(AUC),2); 
            AUC_session{n} = AUC_sorted; % AUC: timeBin x nBoot
        end
        AUC_recSite{session,4} = AUC_session;
        clear AUC_session
    end
    AUC_vSP = [AUC_vSP; AUC_recSite];
    clear AUC_recSite
end
aucPaths = fullfile(mainDir,['vSP_' num2str(nBoot) 'nBoot']);
save(aucPaths, 'AUC_vSP');

%% histogram
% distribution of mean AUC of touch-evoked responses
% label significant units
% normalization: probability
bin = 0.01; % sec
analysisWindow = [0 0.25];
binNum = (analysisWindow(2) - analysisWindow(1))/bin;
% timeWindow = [-1 2.5];
% analysisWindow = [0 0.25];
% binStart =(analysisWindow(1)-timeWindow(1))/bin+1;
% binEnd = (analysisWindow(2)-timeWindow(1))/bin;
save_path= ('E:\CM_NeuralActivity_Analysis\ROC\10ms_bin\decision_probability');
nBoot =1000;
% load('E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity\Unit_AUC_tCorrect_1000nBoot')
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
alpha = 0.05;

for b=1:binNum
    figure('Position',[0 0 800 200])
    for i=1: length(recSites)
        AUC_recSite = AUC_vDP(strcmp(AUC_vDP(:,3), recSites{i}), 4);
        AUC_sessions = [];
        % concatenation across sessions
        for session=1:size(AUC_recSite,1)
            AUC_session = AUC_recSite{session};
            AUC_sessions = [AUC_sessions AUC_session];
        end
        N = size(AUC_sessions,2);
        alpha_Bonferroni = alpha/N; % Bonferroni correction
        AUC_mean = cellfun(@(x) mean(x(b,:)), AUC_sessions,'UniformOutput', false);
        AUC_upCI = cellfun(@(x) x(b,ceil((1-alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        AUC_lowCI = cellfun(@(x) x(b,ceil((alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
        sig_percent = round(sum(isSignificant)/length(isSignificant)*100,1);
        subplot(1,4,i)
        h_all = histogram(cell2mat(AUC_mean),[0:0.05:1], 'FaceColor','white');hold on;
        h_sig = histogram(cell2mat(AUC_mean(isSignificant)),[0:0.05:1], 'FaceColor','black', 'FaceAlpha',0.6);hold off; 
        h_all_normalization = h_all.Values/length(AUC_mean);
        h_sig_normalization = h_sig.Values/length(AUC_mean);
        h_all_prb = histogram('BinEdges',[0:0.05:1], 'BinCounts',h_all_normalization,'FaceColor','white');hold on;
        h_sig_prb = histogram('BinEdges',[0:0.05:1], 'BinCounts',h_sig_normalization,'FaceColor','black', 'FaceAlpha',0.6);hold off; 
        max_count = max([h_all_prb.Values h_sig_prb.Values]);
        min_count = min([h_all_prb.Values h_sig_prb.Values]);
        ylim_max = round(max_count+ 0.1*(max_count - min_count),2);
    %     ylim([0 ylim_max]);
        ylim([0 0.5]);
        xlabel('AUC')
        xticks([0:0.1:1])
        xticklabels({0 [] [] [] [] 0.5 [] [] [] [] 1})
    %     yticks([0:0.05:ylim_max])
        yticks([0:0.1:0.6])
        set(gca, 'box','off','TickDir','out')
        text(0.75,0.3,[num2str(sig_percent), '%'])
        if i==1
            ylabel('Proportion of units')
        end
        title(recSites{i})
        
%         block_preference{i,1} = sum(h_sig_normalization(11:end));% tBlock 
%         block_preference{i,2} = sum(h_sig_normalization(1:10));% vBlock
    end
    startTime = round(analysisWindow(1) + (b-1)*bin,2); 
    endTime =  round(analysisWindow(1) + b*bin,2);
    sgtitle([num2str(startTime) ' - ' num2str(endTime) ' s']);

    % Save fig
    save_path= ('E:\CM_NeuralActivity_Analysis\ROC\10ms_bin\decision_probability\vDP_histogram_10msBin_1000nBoot_Bonferroni'); 
    histFigPath = fullfile(save_path,['vDP_histogram_' num2str(startTime) 'to' num2str(endTime) 's.pdf']);
    print(histFigPath,'-dpdf','-painters','-loose');
end

%% decision probability: lick vs nolick

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% Import setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\ROC\10ms_bin\decision_probability';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
nBoot = 1000;
rng(15)

bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow = [0 0.25];
analysisStart =(analysisWindow(1)-timeWindow(1))/bin+1;
analysisEnd = (analysisWindow(2)-timeWindow(1))/bin;


% DP for tactile trials
% pos_tt = {[1 1 1] [1 1 2] [2 1 1] [2 1 2]}; % lick
% neg_tt = {[1 1 3] [2 1 3]}; % nolick
% DP for visual trials
pos_tt = {[1 2 1] [1 2 2] [2 2 1] [2 2 2]}; % lick
neg_tt = {[1 2 3] [2 2 3]}; % nolick

AUC_vDP = [];
for i = 1: length(recSites)
    u = data_all(strcmp(data_all.recSite, recSites{i}), :);
    AUC_recSite = u{:,1:3};
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1};
        N = size(firingRates_session,1); % number of neurons
        trialNum_session = u{session,5}{1,1};
        parfor n = 1:N % parallel for loop 
            firingRates_pos = [];
            firingRates_neg = [];
            for tt=1:length(pos_tt)
                firingRates_pos_tt = squeeze(firingRates_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3),analysisStart:analysisEnd,...
                    1:trialNum_session(n,pos_tt{tt}(1),pos_tt{tt}(2),pos_tt{tt}(3))));
                firingRates_pos = [firingRates_pos firingRates_pos_tt];
            end  
            for tt=1:length(neg_tt)
                firingRates_neg_tt = squeeze(firingRates_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3),analysisStart:analysisEnd,...
                    1:trialNum_session(n,neg_tt{tt}(1),neg_tt{tt}(2),neg_tt{tt}(3))));                
                firingRates_neg = [firingRates_neg firingRates_neg_tt];
            end
           
            posTrials = num2cell(firingRates_pos,2);
            negTrials = num2cell(firingRates_neg,2);
            posTrialNum = size(firingRates_pos,2);
            negTrialNum = size(firingRates_neg,2);
            binNum = size(posTrials,1);
            % bootstrapped posTrials and negTrials
            [~, posTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), posTrials, 'UniformOutput', false); % Get bootstrapped indexes 
            posTrials_bootstrapped = cellfun(@(x,y) x(y), posTrials, posTrials_bootstrapped_ind, 'UniformOutput', false); % Get bootstrapped values
            posTrials_bootstrapped_matrix = mat2cell(cell2mat(posTrials_bootstrapped),[repmat(posTrialNum,[binNum 1])], [ones([nBoot,1])]); % reshape posTrials_bootstrapped for AUC calculation
            [~, negTrials_bootstrapped_ind] = cellfun(@(x) bootstrp(nBoot,[],x), negTrials, 'UniformOutput', false);  
            negTrials_bootstrapped = cellfun(@(x,y) x(y), negTrials, negTrials_bootstrapped_ind, 'UniformOutput', false);
            negTrials_bootstrapped_matrix = mat2cell(cell2mat(negTrials_bootstrapped),[repmat(negTrialNum,[binNum 1])], [ones([nBoot,1])]); 

            Labels = [ones(posTrialNum,1); zeros(negTrialNum,1)];
            [~,~,~,AUC] = cellfun(@(x,y) perfcurve(Labels,[x(:); y(:)],1),...
                               posTrials_bootstrapped_matrix, negTrials_bootstrapped_matrix, 'uni', 0);
                           
            AUC_sorted = sort(cell2mat(AUC),2); 
            AUC_session{n} = AUC_sorted; % AUC: timeBin x nBoot
        end
        AUC_recSite{session,4} = AUC_session;
        clear AUC_session
    end
    AUC_vDP = [AUC_vDP; AUC_recSite];
    clear AUC_recSite
end
aucPaths = fullfile(mainDir,['vDP_' num2str(nBoot) 'nBoot']);
save(aucPaths, 'AUC_vDP');

%% Heatmap for SP and DP across time

bin = 0.01; % sec
analysisWindow = [0 0.25];
binNum = (analysisWindow(2) - analysisWindow(1))/bin;
nBoot =1000;
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
xTickLables = ['0'; repmat(' ',9,1); '1'; repmat(' ',9,1); '2'; repmat(' ',4,1)]; 
alpha = 0.05;

for b=1:binNum
    for i=1: length(recSites)
        AUC_recSite = AUC_vSP(strcmp(AUC_vSP(:,3), recSites{i}), 4);
        AUC_sessions = [];
        % concatenation across sessions
        for session=1:size(AUC_recSite,1)
            AUC_session = AUC_recSite{session};
            AUC_sessions = [AUC_sessions AUC_session];
        end
        N = size(AUC_sessions,2);
        alpha_Bonferroni = alpha/N; % Bonferroni correction
        AUC_mean = cellfun(@(x) mean(x(b,:)), AUC_sessions,'UniformOutput', false);
        AUC_upCI = cellfun(@(x) x(b,ceil((1-alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        AUC_lowCI = cellfun(@(x) x(b,ceil((alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
        isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
        sig_percent = sum(isSignificant)/length(isSignificant);
        probability_heatmap(i,b) = sig_percent;      
    end
end

figure('Position', [0 0 600 300])
h= heatmap(probability_heatmap, 'Colormap',bone, 'ColorLimits',[0 0.5],...
                'ColorbarVisible', 'on','Title','vSP');
h.XDisplayLabels = xTickLables; h.YDisplayLabels = recSite_names;
h.XLabel = 'Time from stimulus onset (s)';

% Save fig
save_path= ('E:\CM_NeuralActivity_Analysis\ROC\10ms_bin\stimulus_probability'); 
histFigPath = fullfile(save_path,['vSP_' num2str(nBoot) 'nBoot.pdf']);
print(histFigPath,'-dpdf','-painters','-loose');


%% Relationship between modulation of touch-eoked activity and task rule encoding in baseline activity 
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
mainDir = 'E:\CM_NeuralActivity_Analysis\ROC';

% touch-evoked activity: tHit vs. tCR 
load('E:\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity\Unit_AUC_tCorrect_1000nBoot');
bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow_tCorrect = [0 0.05];
binStart =(analysisWindow_tCorrect(1)-timeWindow(1))/bin+1;
binEnd = (analysisWindow_tCorrect(2)-timeWindow(1))/bin;
alpha=0.05;
nBoot=1000;
for i=1: length(recSites)
    AUC_recSite = AUC_tCorrect(strcmp(AUC_tCorrect(:,3), recSites{i}), 4);
    AUC_sessions = [];
    % concatenation across sessions
    for session=1:size(AUC_recSite,1)
        AUC_session = AUC_recSite{session};
        AUC_sessions = [AUC_sessions AUC_session];
    end
    N = size(AUC_sessions,2);
    alpha_Bonferroni = alpha/N; % Bonferroni correction
    AUC_mean = cellfun(@(x) mean(x(binStart:binEnd,:),'all'), AUC_sessions,'UniformOutput', false);
    AUC_upCI = cellfun(@(x) mean(x(binStart:binEnd,ceil((1-alpha_Bonferroni/2)*nBoot))), AUC_sessions,'UniformOutput', false);
    AUC_lowCI = cellfun(@(x) mean(x(binStart:binEnd,ceil((alpha_Bonferroni/2)*nBoot))), AUC_sessions,'UniformOutput', false);
    isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
    AUC_tCorrect_summary{i,1} = recSites{i};
    AUC_tCorrect_summary{i,2} = AUC_mean;
    AUC_tCorrect_summary{i,3} = AUC_upCI;
    AUC_tCorrect_summary{i,4} = AUC_lowCI;
    AUC_tCorrect_summary{i,5} = isSignificant;
end  

clearvars -except recSites recSite_names AUC_tCorrect_summary analysisWindow_tCorrect

% baseline activity: tBlock vs. vBlock
load('E:\CM_NeuralActivity_Analysis\ROC\block_selectivity\Unit_AUC_BP_100msBin_1000nBoot_-1to2p5');
bin = 0.1; % sec
timeWindow = [-1 2.5];
analysisWindow_BP = [-0.1 0];
binStart =round((analysisWindow_BP(1)-timeWindow(1))/bin+1);
alpha=0.05;
nBoot =1000;
for i=1: length(recSites)
    AUC_recSite = AUC_rule(strcmp(AUC_rule(:,3), recSites{i}), 4);
    AUC_sessions = [];
    % concatenation across sessions
    for session=1:size(AUC_recSite,1)
        AUC_session = AUC_recSite{session};
        AUC_sessions = [AUC_sessions AUC_session];
    end
    N = size(AUC_sessions,2);
    alpha_Bonferroni = alpha/N; % Bonferroni correction
    AUC_mean = cellfun(@(x) mean(x(binStart,:)), AUC_sessions,'UniformOutput', false);
    AUC_upCI = cellfun(@(x) x(binStart,ceil((1-alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
    AUC_lowCI = cellfun(@(x) x(binStart,ceil((alpha_Bonferroni/2)*nBoot)), AUC_sessions,'UniformOutput', false);
    isSignificant = cell2mat(cellfun(@(x,y) (0.5-x)*(0.5-y)>0, AUC_upCI, AUC_lowCI,'UniformOutput', false));
    AUC_rule_summary{i,1} = recSites{i};
    AUC_rule_summary{i,2} = AUC_mean;
    AUC_rule_summary{i,3} = AUC_upCI;
    AUC_rule_summary{i,4} = AUC_lowCI;
    AUC_rule_summary{i,5} = isSignificant;
end


figure('Position',[0 0 800 150])
% Correlation
for i=1: length(recSites)
    subplot(1,4,i)
    x = cell2mat(AUC_rule_summary{i,2})';
    y = cell2mat(AUC_tCorrect_summary{i,2})';
%     scatter(x,y, 'b');hold on;
%     [rho,pval] = corr(x,y);
    isSig_tCorrect = AUC_tCorrect_summary{i,5}';
    isSig_rule = AUC_rule_summary{i,5}';
    N = sum(isSig_rule);
    scatter(x(isSig_rule),y(isSig_rule), 20, 'b');hold on;
    [rho,pval] = corr(x(isSig_rule),y(isSig_rule));
    
%     scatter(x(~isSig_rule),y(~isSig_rule), 'b');hold on;
%     [rho,pval] = corr(x(~isSig_rule),y(~isSig_rule));
    
%     scatter(x(isSig_tCorrect),y(isSig_tCorrect), 'b');hold on;
%     [rho,pval] = corr(x(isSig_tCorrect),y(isSig_tCorrect));
    
%     scatter(x(isSig_rule & isSig_tCorrect),y(isSig_rule & isSig_tCorrect), 'b', 'filled');
    
    xlim([0,1]);
    ylim([0 1]);
    xlabel({'Discriminability between block types';...
        ['before stimulus onset (' num2str(analysisWindow_BP(1)) ' to ' num2str(analysisWindow_BP(2)) ')']}, 'FontSize', 6);
    ylabel({'Discriminability between tHit and tCR';...
        ['after stimulus onset (' num2str(analysisWindow_tCorrect(1)) ' to ' num2str(analysisWindow_tCorrect(2)) ')']}, 'FontSize', 6);
    text(0.7,0.2,['r = ' num2str(round(rho,2))], 'FontSize', 6);
    text(0.7,0.3,['n = ' num2str(N)], 'FontSize', 6);
    if pval < 0.001
        text(0.7,0.1,'p < 0.001', 'FontSize', 6);
    else
        text(0.7,0.1,['p = ' num2str(round(pval,3))], 'FontSize', 6);
    end
    xline(0.5, '--k');
    yline(0.5, '--k');
    set(gca, 'box','off','TickDir','out')
    title(recSite_names{i});
end

% 
% save
save_path= ('E:\CM_NeuralActivity_Analysis\ROC\BP_and_tCorrect'); 
% FigPath = fullfile(save_path, ['BP_and_tCorrect_' num2str(analysisWindow_tCorrect(2)*1000) 'ms_all_units']);
% FigPath = fullfile(save_path, ['BP_and_tCorrect_' num2str(analysisWindow_tCorrect(2)*1000) 'ms_non-sig_BP_units']);
FigPath = fullfile(save_path, ['BP_and_tCorrect_' num2str(analysisWindow_tCorrect(2)*1000) 'ms_sig_BP_units_v5']);
% FigPath = fullfile(save_path, ['BP_and_tCorrect_' num2str(analysisWindow_tCorrect(2)*1000) 'ms_sig_tCorrect_units_v2']);
print(FigPath,'-dpdf','-painters','-loose');

% 
% % Pool all brain regions
% figure('Position',[0 0 600 600])
% colors = {[0 1 1] [0 0.5 1] [0 0 1] [0 0 0.5]};
% x = [];
% y = [];
% for i=1: length(recSites)
%     a = cell2mat(AUC_rule_summary{i,2})';
%     b = cell2mat(AUC_tCorrect_summary{i,2})';
%     scatter(a,b,'MarkerFaceColor',colors{i},'MarkerEdgeColor','k');hold on;
%     x = [x;a];
%     y = [y;b];
% end
% 
% [rho,pval] = corr(x,y);
% xlim([0,1])
% ylim([0 1])
% text(0.7,0.2,['r = ' num2str(round(rho,2))])
% text(0.7,0.1,['p = ' num2str(round(pval,4))])
% 
% % Number of neurons showing significant discriminability for (1) tHit vs tCR after stimulus onset, (2) tBlock vs vBlock before stimulus onset 
% % conjunction
% for i=1: length(recSites)
%     x = AUC_rule_summary{i,5};
%     y = AUC_tCorrect_summary{i,5};
%     conj = x+y == 2;
%     sum(conj)
% end
