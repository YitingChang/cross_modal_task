%% recording summary 
% zScore, tHit vs. tCR: isResponsive 
recSite_total = unique(zScore_all.recSite);
for i = 1:length(recSite_total)
    tb = zScore_all(strcmp(zScore_all.recSite, recSite_total{i}), :); % select recording sites
    mouse_number = length(unique(tb.MouseName));
    session_number = length(unique(tb.seshDate));
    responsive_unit_number_sum =0;
    total_unit_number_sum = 0;
    for j = 1:size(tb,1) % loop each session
        responsive_unit_number = sum(tb{j,16}{1,1}); % number of responsive units, the 16th column is 'isResponsive'.
        responsive_unit_number_sum = responsive_unit_number_sum + responsive_unit_number;
        total_unit_number = length(tb{j,16}{1,1});
        total_unit_number_sum = total_unit_number_sum + total_unit_number;
    end
    recording_summary{i,1} = recSite_total{i}; % recording site
    recording_summary{i,2} = mouse_number; % how many mice for that recording site
    recording_summary{i,3} = session_number; % how many sessions for that recording site
    recording_summary{i,4} = responsive_unit_number_sum; % how many stable units
    recording_summary{i,5} = total_unit_number_sum; % how many recorded units 
    
    if i == length(recSite_total) 
        recording_summary{i+1,1} = 'All'; % for all recording sites
        recording_summary{i+1,2} = length(unique(zScore_all.MouseName)); % how many mice for recordings
        recording_summary{i+1,3} = sum(cell2mat(recording_summary(:,3))); % how many sessions in total
        recording_summary{i+1,4} = sum(cell2mat(recording_summary(:,4))); % how many responsive units in total
        recording_summary{i+1,5} = sum(cell2mat(recording_summary(:,5))); % how many recorded units in total 
    end    
end
% Convert to a table
VarNames = {'RecSite', 'MouseNumber', 'SessionNumber', 'ResposiveUnitNumber', 'TotalUnitNumber'};
Recording_summary = cell2table(recording_summary,'VariableNames', VarNames);

writetable(Recording_summary, 'E:\CM_NeuralActivity_Analysis\Recording_summary\Recording_summary_zScore_SiProbe.csv');
%% recording summary
% btInd: > 4, non-drifting
 
recSite_total = unique(btInd_all.recSite);
for i = 1:length(recSite_total)
    tb = btInd_all(strcmp(btInd_all.recSite, recSite_total{i}), :); % select recording sites
    mouse_number = length(unique(tb.MouseName));
    session_number = length(unique(tb.seshDate));
    stable_unit_number_sum =0;
    total_unit_number_sum = 0;
    for j = 1:size(tb,1) % loop each session
        stable_unit_number = length(cell2mat(tb{j,5}{1,1})) - sum(cell2mat(tb{j,5}{1,1})); % number of stable units
        stable_unit_number_sum = stable_unit_number_sum + stable_unit_number;
        total_unit_number = length(cell2mat(tb{j,5}{1,1}));
        total_unit_number_sum = total_unit_number_sum + total_unit_number;
    end
    recording_summary{i,1} = recSite_total{i}; % recording site
    recording_summary{i,2} = mouse_number; % how many mice for that recording site
    recording_summary{i,3} = session_number; % how many sessions for that recording site
    recording_summary{i,4} = stable_unit_number_sum; % how many stable units
    recording_summary{i,5} = total_unit_number_sum; % how many recorded units 
    
    if i == length(recSite_total) 
        recording_summary{i+1,1} = 'All'; % for all recording sites
        recording_summary{i+1,2} = length(unique(btInd_all.MouseName)); % how many mice for recordings
        recording_summary{i+1,3} = sum(cell2mat(recording_summary(:,3))); % how many sessions in total
        recording_summary{i+1,4} = sum(cell2mat(recording_summary(:,4))); % how many stable units in total
        recording_summary{i+1,5} = sum(cell2mat(recording_summary(:,5))); % how many recorded units in total 
    end
end
% Convert to a table
VarNames = {'RecSite', 'MouseNumber', 'SessionNumber', 'StableUnitNumber', 'TotalUnitNumber'};
Recording_summary = cell2table(recording_summary,'VariableNames', VarNames);


writetable(Recording_summary, 'E:\CM_NeuralActivity_Analysis\Recording_summary\Recording_summary_btInd_tetrode.csv');
%% recording summary
% ROC: > 4, non-drifting, isResponsive
 
recSite_total = unique(AUC_all.recSite);
for i = 1:length(recSite_total)
    tb = AUC_all(strcmp(AUC_all.recSite, recSite_total{i}) & strcmp(AUC_all.probability_type, 'stimulus') & ...
        strcmp(AUC_all.trial_type, 'tactile'), :); 
    % select recording sites and get unit information from tactile stimulus probability
    mouse_number = length(unique(tb.MouseName));
    session_number = length(unique(tb.seshDate));
    unit_number_sum =0;
    for j = 1:size(tb,1) % loop each session
        unit_number = size(tb{j,6}{1,1},2);
        unit_number_sum = unit_number_sum + unit_number;
    end
    recording_summary{i,1} = recSite_total{i}; % recording site
    recording_summary{i,2} = mouse_number; % how many mice for that recording site
    recording_summary{i,3} = session_number; % how many sessions for that recording site
    recording_summary{i,4} = unit_number_sum; % how many analyzed units    
    if i == length(recSite_total) 
        recording_summary{i+1,1} = 'All'; % for all recording sites
        recording_summary{i+1,2} = length(unique(AUC_all.MouseName)); % how many mice for recordings
        recording_summary{i+1,3} = sum(cell2mat(recording_summary(:,3))); % how many sessions in total
        recording_summary{i+1,4} = sum(cell2mat(recording_summary(:,4))); % how many analyzed units in total
    end
end
% Convert to a table
VarNames = {'RecSite', 'MouseNumber', 'SessionNumber', 'AnalyzedUnitNumber'};
Recording_summary = cell2table(recording_summary,'VariableNames', VarNames);

writetable(Recording_summary, 'E:\CM_NeuralActivity_Analysis\Recording_summary\Recording_summary_AUC_SiProbe.csv');

%% recording summary
% dPCA and PCA for all trial types (pooled FAs)
% load('E:\CM_NeuralActivity_Analysis\dPCA\dPCA_filter_context\FiringRates_Kmin2_Smooth_NoLick_-1-p15\firingRates_array_SiProbe.mat')

% PCA for all trial types except for TTV & VVT (different FAs)
load('E:\CM_NeuralActivity_Analysis\dPCA\dPCA_filter_context\FiringRates_Kmin2_Smooth_NoLick_-1-p15\FiringRates_no_TTV&VVT\firingRates_array_SiProbe_trialNum')

% number of neurons, sessions, mice
recSites = {'left S1' 'left S2','left wM2' 'left AMM' 'left ALM' 'right ALM' 'left Prt' 'right V1'};
recording_summary = [];
for i = 1: length(recSites)
    u = firingRates_all(strcmp(firingRates_all.recSite, recSites{i}), :);
    num_of_mice = length(unique(u.MouseName));
    num_of_sessions = length(unique(u.seshDate));
    num_of_neurons = 0;
    for session=1:size(u,1)
        firingRate_session = u{session,4}{1,1}; 
        N = size(firingRate_session,1); % number of neurons
        num_of_neurons = num_of_neurons + N;    
    end
    recording_summary_recSite = [{recSites{i}} {num_of_mice} {num_of_sessions} {num_of_neurons}];
    recording_summary = [recording_summary; recording_summary_recSite];
end

% Convert to a table
VarNames = {'recSite' 'num_of_mice' 'num_of_sessions' 'num_of_neurons'};
recording_summary = cell2table(recording_summary,'VariableNames', VarNames);

writetable(recording_summary, 'E:\CM_NeuralActivity_Analysis\Recording_summary\Recording_summary_filter_SiProbe_FAs.csv');

%% recording summary
% number of trials of each trial types in each neuron 

recDir = 'E:\CM_NeuralActivity_Analysis\Recording_summary\Figures_trialNum_distribution';

% PCA for all trial types except for TTV & VVT (different FAs)
load('E:\CM_NeuralActivity_Analysis\dPCA\dPCA_filter_context\FiringRates_Kmin2_Smooth_NoLick_-1-p15\FiringRates_no_TTV&VVT\firingRates_array_SiProbe_trialNum')

% number of neurons, sessions, mice
recSites = {'left S1' 'left S2','left wM2' 'left AMM' 'left ALM' 'right ALM' 'left Prt' 'right V1'};
% concatenate sessions
% Import setting                
S = 2; % number of stimuli: tactile and visual stimuli
D = 3; % number of decisions: lick (including right and left lick) and no lick
C = 2; % number of contexts: respond-to-touch and respond-to-light blocks

recording_summary = [];
for i = 1: length(recSites)
    u = firingRates_all(strcmp(firingRates_all.recSite, recSites{i}), :);
    num_of_mice = length(unique(u.MouseName));
    num_of_sessions = length(unique(u.seshDate));
    maxTrialNum = max(cell2mat(u{:,5}),[],'all');
    firingRates_recSite = [];
    trialNum_recSite = [];
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1};
        trialNum_session = u{session,5}{1,1};
        N = size(firingRates_session,1); % number of neurons
        for n = 1:N
            for s = 1:S
                for d = 1:D
                    for c = 1:C
                        firingRates_session(n,s,d,c,:,trialNum_session(n,s,d,c)+1:maxTrialNum) = nan;
                    end
                end
            end
        end         
        firingRates_recSite = [firingRates_recSite;firingRates_session];
        trialNum_recSite = [trialNum_recSite;trialNum_session];
    end
    num_of_neurons = size(firingRates_recSite,1);
    for s = 1:S
        figure('Position', [0,0, 800, 450]); % set figure size
        for d = 1:D
            for c = 1:C
                trialNum_dist = trialNum_recSite(:,s,d,c);
                subplot(3,2,((d-1)*2+c))
                histogram(trialNum_dist,[0:2:80]);
                xlabel('Number of trials');
                ylabel('Number of units');
                title(['(' num2str(s) ',' num2str(d) ',' num2str(c) ')']);
            end
        end
        switch s
            case 1
                titleName = [recSites{i} ' tStim'];
                sgtitle(titleName)
            case 2
                titleName = [recSites{i} ' vStim'];
                sgtitle(titleName)
        end
        % Save figure  
        FigPath = fullfile(recDir,[titleName '.pdf']);
        print(FigPath,'-dpdf','-painters','-loose');
    end
end    
       
%% recording summary
% number of trials of early and late VTTs in each neuron 

recDir = 'E:\CM_NeuralActivity_Analysis\Recording_summary\Figures_trialNum_distribution';

% Map early and late VTTs onto tHit space 
load('E:\CM_NeuralActivity_Analysis\dPCA\dPCA_filter_context\FiringRates_Kmin2_Smooth_NoLick_-1-p15\FiringRates_no_TTV&VVT\firingRates_array_SiProbe_trialNum')

% number of neurons, sessions, mice
recSites = {'left S1' 'left S2','left wM2' 'left AMM' 'left ALM' 'right ALM' 'left Prt' 'right V1'};
% concatenate sessions
% Import setting                
S = 2; % number of stimuli: tactile and visual stimuli
D = 3; % number of decisions: lick (including right and left lick) and no lick
C = 2; % number of contexts: respond-to-touch and respond-to-light blocks

recording_summary = [];
for i = 1: length(recSites)
    u = firingRates_all(strcmp(firingRates_all.recSite, recSites{i}), :);
    MouseNames = [];
    seshDates = [];
    num_of_neurons = 0;
    Early_VTT = [];
    Late_VTT = [];
    for session=1:size(u,1)
        firingRates_session = u{session,4}{1,1};
        trialNum_vBlock = u{session,7}{1,1};
        N = size(firingRates_session,1); % number of neurons                
        
        [~,~,trialNum_VTT] = find(trialNum_vBlock(1,1,1,2,:));
        isEarly_VTT = trialNum_VTT<=10;
        isLate_VTT = trialNum_VTT>10;
        VTTs_session = [isEarly_VTT isLate_VTT];
        
        if sum(sum(VTTs_session,1)>=1)==2 
           MouseName = u{session,1}{1,1}; 
           MouseNames = [MouseNames; {MouseName}]; 
           seshDate = u{session,2}{1,1}; 
           seshDates = [seshDates; {seshDate}];
           num_of_neurons = num_of_neurons + N; 
           Early_VTT_session = repmat(sum(isEarly_VTT),[N 1]); 
           Early_VTT = [Early_VTT; Early_VTT_session];
           Late_VTT_session = repmat(sum(isLate_VTT),[N 1]);
           Late_VTT = [Late_VTT; Late_VTT_session];
        end
    end
    
    num_of_mice = length(unique(MouseNames));
    num_of_sessions = length(unique(seshDates));
    recording_summary_recSite = [{recSites{i}} {num_of_mice} {num_of_sessions} {num_of_neurons}];
    recording_summary = [recording_summary; recording_summary_recSite];
    
    figure('Position', [0,0, 750, 300]); 
    subplot(1,2,1)
    histogram(Early_VTT, [0:1:15]);
    xlabel('Number of trials');
    ylabel('Number of neurons');
    title('early')
    subplot(1,2,2)
    histogram(Late_VTT, [0:1:15]);
    xlabel('Number of trials');
    ylabel('Number of neurons');
    title('late')
    sgtitle([recSites{i} ' tFA(rule error)'])
    
    % Save figure  
    FigPath = fullfile(recDir,[recSites{i} '_VTT_early&late' '.pdf']);
    print(FigPath,'-dpdf','-painters','-loose');
    
end  

% Convert to a table
VarNames = {'recSite' 'num_of_mice' 'num_of_sessions' 'num_of_neurons'};
recording_summary = cell2table(recording_summary,'VariableNames', VarNames);

writetable(recording_summary, 'E:\CM_NeuralActivity_Analysis\Recording_summary\Recording_summary_filter_SiProbe_VTT_early&late.csv');

