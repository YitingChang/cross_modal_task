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

%% Determine tHit-tCR selectivity when three consecutive time bins are significant. 
mainDir = 'E:\Cross-Modal_Project\Revision\tHit-tCR_selectivity_and_subspace_overlap';
bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow = [0 0.15];
binStart =round((analysisWindow(1)-timeWindow(1))/bin)+1;
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);
nBoot =1000;
load('E:\Cross-Modal_Project\CM_NeuralActivity_Analysis\ROC\tHit_tCR_selectivity\Unit_AUC_tCorrect_1000nBoot')
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
alpha = 0.05;
significant_bins = 3;
% 
AUC_tCorrect_summary =[];    
for i=1: length(recSites)
    auc_recSite = AUC_tCorrect(strcmp(AUC_tCorrect(:,3), recSites{i}), 4);
    session_info_recSite = AUC_tCorrect(strcmp(AUC_tCorrect(:,3), recSites{i}), 1:3);
    AUC_sessions = [];
    session_names = [];
    % concatenation across sessions
    for session=1:size(auc_recSite,1)
        auc_session = auc_recSite{session};
        AUC_sessions = [AUC_sessions auc_session];
        num_neurons = size(auc_session,2);
        session_name = repmat(session_info_recSite(session,2),1,num_neurons);
        session_names = [session_names session_name];
    end
    N = size(AUC_sessions,2);
    alpha_Bonferroni = alpha/N; % Bonferroni correction
    
    AUC_mean = cellfun(@(x) mean(x(binStart:binEnd,:),'all'), AUC_sessions,'UniformOutput', false);
    for unit=1:N
        AUC_unit = AUC_sessions{unit};
        for b= binStart:binEnd
            AUC_upCI = AUC_unit(b,ceil((1-alpha_Bonferroni/2)*nBoot));
            AUC_lowCI = AUC_unit(b,ceil((alpha_Bonferroni/2)*nBoot));
            isSignificant_bin(b-binStart+1) = (0.5-AUC_upCI)*(0.5-AUC_lowCI)>0;
        end
        % If three consecutive time bins are significant, this unit show significant selectivity.
        isSignificant_bin_matrix = [isSignificant_bin; [isSignificant_bin(2:end), 0]; [isSignificant_bin(3:end), 0, 0]];
        isSignificant_three_bins = ismember(sum(isSignificant_bin_matrix,1),significant_bins);
        isSignificant(unit) = sum(isSignificant_three_bins)>0 ;
    end
    
    for session=1:size(auc_recSite,1)
        session_name = session_info_recSite{session,2};
        % Select neurons from the current session
        isSession = strcmp(session_names, session_name);
        session_info_recSite{session, 4} = sum(isSession); % number of neurons in the current session
        session_info_recSite{session, 5} = sum(isSession & isSignificant);  % number of siginificant neurons in the current session
        session_info_recSite{session, 6} = sum(isSession & isSignificant)/sum(isSession); % percentage of neurons in the current session
        session_info_recSite{session, 7} = AUC_mean(isSession); % AUC mean for each neurons in the current session
        session_info_recSite{session, 8} = isSignificant(isSession);
    end  
    AUC_tCorrect_summary = [AUC_tCorrect_summary; session_info_recSite];  
    clearvars isSignificant isSignificant_bin
end
% Save
AUC_tCorrect_summary = cell2table(AUC_tCorrect_summary,'VariableNames',{'mouse_name', 'session_name', 'recSite',...
                                                        'num_neurons', 'num_sig_neurons', 'perc_sig_neurons','AUC_mean', 'isSignificant'});
save_path = fullfile(mainDir,'Unit_AUC_tCorrect_summary_3bins');
save(save_path, 'AUC_tCorrect_summary');

%% histogram
% distribution of mean AUC of touch-evoked responses
% label significant units
% normalization: probability
mainDir = 'E:\Cross-Modal_Project\Revision\tHit-tCR_selectivity_and_subspace_overlap';
bin = 0.01; % sec
timeWindow = [-1 2.5];
analysisWindow = [0 0.15];
binStart =round((analysisWindow(1)-timeWindow(1))/bin)+1;
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);
figure('Position',[0 0 200 800])
nBoot =1000;
load(fullfile(mainDir,'Unit_AUC_tCorrect_summary_3bins'))
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
alpha = 0.05;

for i=1: length(recSites)
    AUC_tCorrect_recSite = AUC_tCorrect_summary(strcmp(AUC_tCorrect_summary.recSite, recSites{i}),:);
    AUC_mean = [];
    isSignificant = [];
    % concatenation across sessions
    for session=1:size(AUC_tCorrect_recSite,1)
        AUC_session = AUC_tCorrect_recSite.AUC_mean{session};
        AUC_mean = [AUC_mean AUC_session];
        isSignificant_session = AUC_tCorrect_recSite.isSignificant{session};
        isSignificant = [isSignificant isSignificant_session];
    end
    sig_percent = round(sum(isSignificant)/length(isSignificant)*100,1); 
    
    subplot(4,1,i)
    h_all = histogram(cell2mat(AUC_mean),[0:0.05:1], 'FaceColor','white');hold on;
    h_sig = histogram(cell2mat(AUC_mean(logical(isSignificant))),[0:0.05:1], 'FaceColor','black', 'FaceAlpha',0.6);hold off; 
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
histFigPath = fullfile(mainDir, ['AUC_histogram_tCorrect_Bonferroni_' num2str(nBoot) 'nBoot_p0top15_3SigBins.pdf']);
print(histFigPath,'-dpdf','-painters','-loose');

