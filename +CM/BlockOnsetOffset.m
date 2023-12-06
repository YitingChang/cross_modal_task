% Get block onset and offset for each session
% Convert block onset/offset to index
% debugging: some MSessionExplore do not have reference time

load('E:\seFilePaths_all');

for m = 21 % 1:size(seFilePaths_all,1)
    tic
    % Get block ind for individual mice
    MouseName = seFilePaths_all.MouseName{m};
    seFilePaths = seFilePaths_all.seFilePaths{m};
    
    BlockTimes = [];
    for i = 1: length(seFilePaths)
        load(seFilePaths{i})
        
        % Get block index
        [~,~,~,~,tBlockInd, vBlockInd] = CM.QualityControl.GetBlockInfo(se);  
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

        % Reference times are stimulus onsets. Stimulus onsets are used as trial onsets. 
        reference_time = se.GetReferenceTime('spikeTime'); 
        tBlock_start_time = reference_time(tBlock_start_trial).';
        tBlock_end_time = reference_time(tBlock_end_trial).';
        vBlock_start_time = reference_time(vBlock_start_trial).';
        vBlock_end_time = reference_time(vBlock_end_trial).';
        
        % Load session information 
        seshDate = se.userData.sessionInfo.seshDate;
        recSite = se.userData.sessionInfo.recSite;
        
        % Concatenation 
        new_BlockTimes = [{MouseName} {seshDate} {recSite}...
                            tBlock_start_time tBlock_end_time vBlock_start_time vBlock_end_time];
%                            {tBlock_start_time} {tBlock_end_time} {vBlock_start_time} {vBlock_end_time}];
        BlockTimes = [BlockTimes; new_BlockTimes];
        
        clearvars -except seFilePaths_all seFilePaths BlockTimes MouseName
    end
    
    % Convert to a table
    VarNames = {'MouseName' 'seshDate' 'recSite'...
            'tBlock_start_time' 'tBlock_end_time' 'vBlock_start_time' 'vBlock_end_time'};
    BlockTimes = cell2table(BlockTimes,'VariableNames', VarNames);
    
    toc
    % Save block onset and offset for each mouse
    BlockTimesPath = fullfile('E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\BlockTimes', ...
                        [MouseName,'_BlockTimes']);
    save(BlockTimesPath, 'BlockTimes');
end

%% Concatenation
% Convert table to cell
% Python cannot read matlab table
BlockTimesPaths = MBrowse.Files('E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\BlockTimes');
BlockTimes_all = [];
for i = 1: length(BlockTimesPaths)
    data = load(BlockTimesPaths{i});
    BlockTimes_mouse = data.BlockTimes;
    if iscell(BlockTimes_mouse)
        BlockTimes_all = [BlockTimes_all; BlockTimes_mouse];
    else
        BlockTimes_mouse = table2cell(BlockTimes_mouse);
        BlockTimes_all = [BlockTimes_all; BlockTimes_mouse];
    end
end
% Save 
save('E:\CM_NeuralActivity_Analysis\Unit_quality_metrics\BlockTimes\BlockTimes_SiProbe', 'BlockTimes_all', '-v7');