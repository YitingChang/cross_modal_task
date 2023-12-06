% load data array
load('E:\Data_array\AllTrials_100bin_noSmooth\data_array.mat')
mainDir = 'E:\CM_Behavior_Analysis\firstHit_distribution';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
B=2; % number of blocks: respond-to-touch and respond-to-light blocks
S=2; % number of stimuli: tactile and visual stimuli
D=3; % number of decisions: lick (including right and left lick) and no lick

%% histogram of trial numbers for the first hit trial (all recSites)    
% tHit vs. vHit
titles = {'Respond-to-touch' 'Respond-to-light'};
u = data_all(ismember(data_all.recSite, recSites), :);
figure('Position',[0 0 750 300])
for tt = 1:2
    subplot(1,2,tt)
    trialNum_firstHit = [];
    for session=1:size(u,1)
        v = u{session,12};
        switch tt
            case 1 % tHit
                trialNum_firstHit_session = v.trialNum_tHit;
            case 2 % vHit
                trialNum_firstHit_session = v.trialNum_vHit;
        end
        trialNum_firstHit = [trialNum_firstHit trialNum_firstHit_session];
    end
%     h = histogram(cell2mat(trialNum_firstHit),[1:2:24],'FaceColor',[0.5 0.5 0.5]); 
%     ylim([0 25]);
%     text(9.5,22,'switch cue');
%     ylabel('Number of blocks');
    h = histogram(cell2mat(trialNum_firstHit),[1:2:24],'Normalization','probability','FaceColor',[0.5 0.5 0.5]); 
    ylim([0 0.25]);
    text(9.5,0.23,'switch cue')
    ylabel('Fraction of blocks');
    xlim([0 24]);
    set(gca, 'box','off', 'TickDir','out');
    xlabel('Trial number of first hit after block switch');
    xline(9,'--k');
    title(titles{tt});
    % fraction of first hit trials occouring before the switch cue
    text(5,0.15,num2str(round(sum(h.Values(1:4)),2)));
end

% Save figure  
FigPath = fullfile(mainDir,'TrialNum_of_firstHit_fraction_without_1stBlock_v2');
print(FigPath,'-dpdf','-painters','-loose');
 