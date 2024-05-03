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

%% Heatmap (tBlock vs. vBlock)
% initial activity of correct trials 
% zScore:(x- mean(baseline))/std(baseline) 
% initial activity: -100 ~ 0 ms
% Chang et al., Fig 2b

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

