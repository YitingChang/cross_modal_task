% Fractions of trial outcomes

% Load data
load('E:Data_array\AllTrials_10bin_Smooth\data_array.mat')
% Setting 
save_path = 'E:\CM_Behavior_Analysis\Behavioral_performance';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
% Brain area selection 
u = data_all(ismember(data_all.recSite, recSites),:);

mouseNames = unique(u.MouseName);
behav_summary = [];
for mouse = 1:length(mouseNames)
    v = u.behavPerformance(strcmp(u.MouseName, mouseNames{mouse}), :);
    if size(v,1)>1
        MEAN_mouse = mean(v{:,:},1);
    else
        MEAN_mouse = v{:,:};
    end
    behav_summary =[behav_summary; MEAN_mouse];
end

MEAN_all_tt = [];
SEM_all_tt = [];
for tt=1:size(behav_summary,2)
    k = behav_summary(:,tt);
    MEAN_all = mean(k);
    SEM_all = std(k)/sqrt(length(k));
    MEAN_all_tt = [MEAN_all_tt; MEAN_all];
    SEM_all_tt = [SEM_all_tt; SEM_all];
end

%% Stacked bar plot for trial outcomes in respond-to-touch and respond-to-light blocks
% Chang et al. Figure 1d

figure('Position', [0,0, 400, 600]); % set figure size
tBlock_mean = MEAN_all_tt(1:4); 
tBlock_sem = SEM_all_tt(1:4);
vBlock_mean = MEAN_all_tt(5:8); 
vBlock_sem = SEM_all_tt(5:8); 
b = bar([tBlock_mean';vBlock_mean'],'stacked', 'FaceColor', 'flat');
bar_colors = {[0 0 1], [1 0 0], [0.8 0.8 0.8], [0.5 0.5 0.5]};
for n = 1:4
    b(:,n).CData = bar_colors{n};
    text(0.75,1.1+0.05*n, [num2str(round(tBlock_mean(n),2)) ' (' num2str(round(tBlock_sem(n),2)) ')']);
    text(1.75,1.1+0.05*n, [num2str(round(vBlock_mean(n),2)) ' (' num2str(round(vBlock_sem(n),2)) ')']);
end
hold on
errorbar(cumsum([tBlock_mean vBlock_mean])', [tBlock_sem'; vBlock_sem'], '.k');
xlim([0.5 3.5]);
ylim([0 1.5]);
yticks(0:0.2:1)
ylabel('Fraction of trials')
xticklabels(['Respond-to-touch'; 'Respond-to-light'])
legend([{'Hit'} {'CR'} {'Miss'} {'FA'}].','location','northeast')
set(gca, 'box','off','TickDir','out')

% Test respond-to-touch and respond-to-light performance
tBlock_perf = behav_summary(:,17); % respond-to-touch performance
vBlock_perf = behav_summary(:,18); % respond-to-light performance
[h,p,ci,stats] = ttest(tBlock_perf,vBlock_perf); % paired-sample t-test

for b=1:2
    text(b-0.25, 1.45, 'correctRate = ')
    text(b-0.25, 1.4, [num2str(round(MEAN_all_tt(16+b),2)) ' (' num2str(round(SEM_all_tt(16+b),2)) ')'])
end
text(2.7,0.6,{'paired ttest',['p-value= ' num2str(round(p,2))]})

%% Save fig
behavFigPath = fullfile(save_path,'Behavioral_performance.pdf');
print(behavFigPath,'-dpdf','-painters','-loose');

%% Stacked bar plot for FA types in respond-to-touch and respond-to-light blocks
FA_tBlock_mean = MEAN_all_tt(9:11)/MEAN_all_tt(4);
FA_vBlock_mean = MEAN_all_tt(12:14)/MEAN_all_tt(8);
FA_tBlock_sem = SEM_all_tt(9:11)/MEAN_all_tt(4);
FA_vBlock_sem = SEM_all_tt(12:14)/MEAN_all_tt(8);

subplot(1,2,2)
b_FA = bar([FA_tBlock_mean, FA_vBlock_mean].','stacked', 'FaceColor', 'flat');
bar_colors = {[0 1 1], [1 0 1], [1 1 0]};
for n = 1:3
    b_FA(:,n).CData = bar_colors{n}; 
end
hold on
errorbar(cumsum([FA_tBlock_mean, FA_vBlock_mean])', [FA_tBlock_sem, FA_vBlock_sem].', '.k');
xlim([0.5 3.5]);
ylim([0 1.1]);
ylabel('Fraction of trials')
xticklabels(['Respond-to-touch'; 'Respond-to-light'])
legend([{'cl'} {'re'} {'ex'}].','location','northeast') % compulsive licking, rule error, exploration

