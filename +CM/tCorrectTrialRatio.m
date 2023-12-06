%% GetAmplitudeRatio(MouseName, seDir)
% Amplitude ratios from spike rate
load('E:\seFilePaths_all');
tCorrectTrialRatio_summary =[];
for m = 1:size(seFilePaths_all,1)
    tic
    MouseName = seFilePaths_all.MouseName{m};
    seFilePaths = seFilePaths_all.seFilePaths{m};

    for i = 1: length(seFilePaths)
        load(seFilePaths{i})
        tCorrectTrialRatio = [];
        spkRs_response = CM.Analysis.GetSpikeRate(se, 0.25, 0, 0.25); % response = 0 ~ 0.25 sec
        trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
             'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
        ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
        'VariableNames', trialTypes);
    
        isResponsive = CM.Analysis.GetResponsiveUnit(se, [0,2]);        
         % select unit
        for u = 2:size(spkRs_response, 2) % first column is time and unit starts from 2nd column
            if isResponsive(u-1) 
                sr = cell2mat(spkRs_response{:,u});
                % select trial types
                sr_TTT = sr(ttInd.TTT);
                sr_VTN = sr(ttInd.VTN);
                % Haldane's correction: adding 0.5 to all cells
                ratio = (mean(sr_TTT)+0.5)/(mean(sr_VTN)+0.5);
                % summary
                tCorrectTrialRatio(u-1) = ratio;
            else
                tCorrectTrialRatio(u-1) = NaN;
            end
        end    

        if ~isempty(tCorrectTrialRatio)
            % Load session information 
            MouseName = se.userData.sessionInfo.MouseName;
            seshDate = se.userData.sessionInfo.seshDate;
            recSite = se.userData.sessionInfo.recSite;

            % 
            new_tCorrectTrialRatio = [{MouseName} {seshDate} {recSite} {tCorrectTrialRatio}];

            % Concatenation 
            tCorrectTrialRatio_summary = [tCorrectTrialRatio_summary; new_tCorrectTrialRatio];
        else
            % Stop if there is no responsive unit in this session
            continue
        end
        clearvars -except seFilePaths_all seFilePaths tCorrectTrialRatio_summary MouseName
    end
    toc
end

% Convert to a table
VarNames = {'MouseName', 'seshDate', 'recSite', 'tCorrectTrialRatio'};
tCorrectTrialRatio_summary = cell2table(tCorrectTrialRatio_summary,'VariableNames', VarNames);

% Save tCorrectTrials_summary summary for each mouse 
tCorrectPaths = fullfile('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Amplitude\tCorrectAmpRatio_SiProbe');
save(tCorrectPaths, 'tCorrectTrialRatio_summary');

%% Ratio of tHit and tCR zScore
recSite = {'left S1', 'left S2', 'left wM2', 'left AMM','left ALM'};
save_path = 'E:\CM_NeuralActivity_Analysis\tCorrectTrials\Amplitude';

for i = 1: length(recSite)
    u = tCorrectAmp_all(strcmp(tCorrectAmp_all.recSite, recSite{i}), :); % select recording sites
    MouseNames = unique(u.MouseName);
    for mouse = 1:length(MouseNames)
        m = u(strcmp(u.MouseName, MouseNames(mouse)),:); 
        tHit = [];
        tCR = [];
        for session = 1:size(m,1) % loop each session
            s = cell2mat(m{session,4}{1,1}); % get tCorrect summary

            % remove unit having NaN and inf
            s(:,any(isnan(s), 1)) = [];
            s(:,any(isinf(s), 1)) = [];
 
            % Get all units
            tHit_all = s(1,:);
            tCR_all = s(2,:);
            tHit = [tHit tHit_all];
            tCR = [tCR tCR_all];
        end
        ratio = mean(tCR)/mean(tHit);
        tCorrectRatio_mouse{i,mouse} = ratio;
    end
    y = cell2mat(tCorrectRatio_mouse(i,:));
    [h, p] = ttest(y,1); % one sample t test
    tCorrectRatio_ttest{i,1} = mean(y);
    tCorrectRatio_ttest{i,2} = std(y)/ sqrt(size(y,1));
    tCorrectRatio_ttest{i,3} = length(y); 
    tCorrectRatio_ttest{i,4} = p;
    
end

% plot
figure('Position', [0,0, 400, 300]); % set figure size
for i = 1: length(recSite)
    y = cell2mat(tCorrectRatio_mouse(i,:));
    x = repmat(i,1,length(y));
    scatter(x,y,10,'k', 'jitter','on', 'jitterAmount',0.1)
    hold on
    errorbar(i-0.3,mean(y), std(y)/ sqrt(size(y,1)),...
    'o','MarkerSize',6,'Color','blue', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
end
xlim([0.5 length(recSite)+0.5]);
ylim([-0.2 1.2]);
ylabel('tCR/tHit')
xticklabels({'wS1' 'wS2' 'wM2' 'AMM' 'tjM2'})
set(gca,'xtick',1:1:8)
% Save fig
tCorrectAmpFigPath = fullfile(save_path,'tCorrectAmpRatio_mixed.pdf');
print(tCorrectAmpFigPath,'-dpdf','-painters','-loose')

%% Amplitude ratios from zScore - 2020/07/10

load('E:\CM_NeuralActivity_Analysis\zScore\zScore_qualityFilter_NoLickWindow\zScore_SiProbe.mat');

% setting
bin = 0.01;
traceWindow = [-1 2.5]; 
analysisWindow = [0 0.25];
timeBin_start = (analysisWindow(1)-traceWindow(1))*(1/bin) + 1;
timeBin_end = (analysisWindow(2)-traceWindow(1))*(1/bin);
save_path = 'E:\CM_NeuralActivity_Analysis\zScore\zScore_qualityFilter_NoLickWindow\tCorrect\';

tCorrectTrialMean_all = [];
for session = 1:size(zScore_all,1)
    MouseName = zScore_all.MouseName{session};
    seshDate = zScore_all.seshDate{session};
    recSite = zScore_all.recSite{session};
    
    tHit = zScore_all.TTT{session};
    tCR =  zScore_all.TTN{session};
    
    tHit_mean = mean(tHit(timeBin_start:timeBin_end,:),1);
    tCR_mean = mean(tCR(timeBin_start:timeBin_end,:),1);
        
    new_tCorrectTrialMean = [{MouseName} {seshDate} {recSite} {tHit_mean} {tCR_mean}];
    
    % Concatenation 
    tCorrectTrialMean_all = [tCorrectTrialMean_all; new_tCorrectTrialMean];
end
% Convert to a table
VarNames = {'MouseName', 'seshDate', 'recSite', 'tHit', 'tCR'};
tCorrectTrialMean_all = cell2table(tCorrectTrialMean_all,'VariableNames', VarNames);

% Plot - all neurons
recSite = {'left S1', 'left S2', 'left wM2', 'left ALM'};
figure('Position', [0 0 800 800])
for i = 1: length(recSite)
    u = tCorrectTrialMean_all(strcmp(tCorrectTrialMean_all.recSite, recSite{i}), :); % select recording sites
    x = cell2mat(u.tHit');
    y = cell2mat(u.tCR');
    subplot(2,2,i)
    scatter(x,y,'b'); hold on
    plot([-10 20], [-10 20], '--k');
    xlim([-2 10]); ylim([-2 10]);
    xlabel('tHit zScore');ylabel('tCR zScore');
    title(recSite{i})
end

% Save figure
tCorrectZScoreFigPath = fullfile(save_path,'tCorrect_zScore.pdf');
print(tCorrectZScoreFigPath,'-dpdf','-painters','-loose')


% Plot - averages for individual mice
figure('Position',[0 0 400 400])
for i = 1: length(recSite)
    u = tCorrectTrialMean_all(strcmp(tCorrectTrialMean_all.recSite, recSite{i}), :); % select recording sites
    MouseNames = unique(u.MouseName);
    ratios = [];
    for mouse = 1:length(MouseNames)
        m = u(strcmp(u.MouseName, MouseNames(mouse)),:); 
        ratio = mean(cell2mat(m.tCR'))/ mean(cell2mat(m.tHit'));
        ratios = [ratios ratio];
    end
    x = repmat(i,1,length(ratios));
    scatter(x,ratios,'k', 'jitter','on', 'jitterAmount',0.1);hold on
    errorbar(i-0.3,mean(ratios), std(ratios)/ sqrt(size(ratios,1)),...
    'o','MarkerSize',6,'Color','blue', 'MarkerEdgeColor','blue','MarkerFaceColor','blue');
end

xlim([0.5 length(recSite)+0.5]);
ylim([0 1]);
ylabel('tCR/tHit')
xticklabels({'wS1' 'wS2' 'wM2' 'tjM2'})
set(gca,'xtick',1:1:8)
% Save fig
tCorrectRatioFigPath = fullfile(save_path,'tCorrect_Ratio.pdf');
print(tCorrectRatioFigPath,'-dpdf','-painters','-loose')

