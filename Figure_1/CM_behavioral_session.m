%%
% top: blcok types (respond-to-touch: blue; respond-to-light: red)
% middle: tactile stimuli 
% bottom: visual stimuli 
% right lick (green); left lick (orange); no lick (gray)
% block start (dash line)
% cue for block transition (star)

% setting
y_idx = [2 1];
lick_color = {[0.5 0.2 0.55] [1 0.4 0.15] [0.75 0.75 0.75]};
save_path = 'E:\CM_Behavior_Analysis\Behaivoral_session\';
% load('E:\seFilePaths_SiProbe'); % File paths for MSessionExplorers
% Directory
Mice = {'YT084'};
Disks = {'F:\'};
m=1;
mouseName = Mice{m};
MainDir = Disks{m};
seDir = [MainDir, mouseName, '\MSessionExplorer_Ephy'];
cd(seDir)
seFileLists = struct2cell(dir('*se*'));
seFileNames = seFileLists(1,:);
seFilePaths = cellfun(@(x) fullfile(seDir, x),seFileNames,'UniformOutput',false);
%%
for i = 1: length(seFilePaths)
    load(seFilePaths{i})
    
    mouseName = se.userData.sessionInfo.MouseName;
    seshDate = se.userData.sessionInfo.seshDate;
    recSite = se.userData.sessionInfo.recSite;
    
    isTblock = ismember(se.GetColumn('behavValue', 'blockType'), 'Whisker');
    isVblock = ismember(se.GetColumn('behavValue', 'blockType'), 'Visual');
    isTstim = ~isnan(cell2mat(se.GetColumn('behavTime', 'somOnset')));
    isVstim = ~isnan(cell2mat(se.GetColumn('behavTime', 'visOnset'))); 
    stim_idx = {isTstim isVstim};
    response = se.GetColumn('behavValue', 'response');
    isRlick = logical(cell2mat(cellfun(@(x) x==1, response, 'UniformOutput',false))); 
    isLlick = logical(cell2mat(cellfun(@(x) x==2, response, 'UniformOutput',false)));
    isNolick = logical(cell2mat(cellfun(@(x) x==0, response, 'UniformOutput',false)));
    lick_idx = {isRlick isLlick isNolick};

    isChange = ismember(se.GetColumn('behavValue', 'blockType'), 'Change') ;
    switch_cue_idx = find(isChange == 1);
    tBlockInd = double(isTblock);
    vBlockInd = double(isVblock);
    % Assign block type to the cahnge trials  
    for t = 2: length(tBlockInd)
        if isChange(t) == 1 && tBlockInd(t-1) ==1
            tBlockInd(t) = 1;
        elseif isChange(t) == 1 && vBlockInd(t-1) ==1
            vBlockInd(t) = 1;
        end 
    end
    % Find first and last trial of blocks
    [tBlock_start_trial,~] = find(diff([0;tBlockInd])== 1); 
    [vBlock_start_trial,~] = find(diff([0;vBlockInd])== 1);
    if tBlockInd(end) == 1
        [tBlock_end_trial,~] = find(diff([tBlockInd;0])== -1); 
        [vBlock_end_trial,~] = find(diff([vBlockInd])== -1);
    else
        [tBlock_end_trial,~] = find(diff([tBlockInd])== -1); 
        [vBlock_end_trial,~] = find(diff([vBlockInd;0])== -1);
    end        
    NumOfTrials = length(isTblock);

    figure('Position',[0 0 800 250])
    % block type 
    x1 = find(tBlockInd == 1);
    y1 = repmat(3.3,[sum(tBlockInd),1]);
    scatter(x1,y1,12,'b','s','filled'); hold on;
    x2 = find(vBlockInd == 1);    
    y2 = repmat(3.3,[sum(vBlockInd),1]);
    scatter(x2,y2,12,'r','s','filled');

    % lick responses for tactile(s=1) and visual(s=2) trials 
    for s=1:2
        stim = stim_idx{s};
        y = y_idx(s);
        for l=1:3
            lick = lick_idx{l};
            is_stim_lick = (stim+lick) ==2;
            x = find(is_stim_lick == 1);
            y_start = repmat(y,[sum(is_stim_lick),1]);
            y_stop = repmat(y+0.7,[sum(is_stim_lick),1]);
            plot([x.';x.'],[y_start.';y_stop.'],'Color', lick_color{l})
        end
        x_cue = switch_cue_idx(stim(switch_cue_idx)==1);
        y_cue = repmat(y+0.8, [length(x_cue),1]);
        scatter(x_cue, y_cue, 16, 'ok', 'filled')
    end

    % block start
    block_start = sort([tBlock_start_trial; vBlock_start_trial]);
    x = block_start(2:end)-0.5; % -0.5 to avoid cover the first trials
    y_start = repmat(min(y_idx)-0.2,[length(x),1]);
    y_stop = repmat(max(y_idx)+1.5,[length(x),1]);
    plot([x.';x.'],[y_start.';y_stop.'],'--k','Linewidth',1)


    set(gca, 'box','off','TickDir','out');
    xlim([0 NumOfTrials])
    ylim([0.8 3.5])
    xticks([]);
    xlabel('Time within a session')
    yticks([1.4 2.4 3.3]);
    yticklabels({'Visual' 'Tactile' 'Rule'})

    title([mouseName '\_' seshDate '\_' recSite])

    % Save fig
    fig = gcf;
    lickFigPath = fullfile(save_path,[mouseName '_' seshDate '_' recSite '.pdf']);
    print(lickFigPath,'-dpdf','-painters','-loose');    
end