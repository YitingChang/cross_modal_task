% Matlab structure to Python dictionary

% all_beh in python
% Index(['block_type', 'cluster_name', 'correct', 'date', 'first_lick',
%        'identified', 'last_lick', 'licks_left', 'licks_right', 'mouse_name',
%        'response', 'session_type', 'spike_times', 'spike_times(stim_aligned)',
%        'stim_duration', 'stim_offset', 'stim_onset', 'touch_stimulus',
%        'trial_num', 'trial_onset', 'trial_type', 'uni_id', 'vis_stimulus',
%        'opto_stim_onsets', 'inhibition_session'],
%       dtype='object')
     
%%
clearvars -except se
% Load behavValue and behavTime
behavValue = se.GetTable('behavValue'); 
behavValue = behavValue(:,1:6);
behavTime = se.GetTable('behavTime');

% Get correct responses
TTT = se.GetColumn('behavValue', 'TTT');
VVV = se.GetColumn('behavValue', 'VVV');
TVN = se.GetColumn('behavValue', 'TVN');
VTN = se.GetColumn('behavValue', 'VTN');

correct = num2cell(TTT | VVV | TVN | VTN); % correct responses

% Add session info to each trials
MouseName = mat2cell(se.userData.sessionInfo.MouseName,1);
seshDate = mat2cell(se.userData.sessionInfo.seshDate,1);
seshType = mat2cell(se.userData.sessionInfo.seshType,1);
mouse_name = repmat(MouseName, length(behavValue.bct_trialNum),1);
sesh_date = repmat(seshDate,length(behavValue.bct_trialNum),1);
sesh_type = repmat(seshType,length(behavValue.bct_trialNum),1);

for i  = 1: length(sesh_type)
    if strcmp(sesh_type{i},'sham')
        inh_session{i,1} = 'sham';
    else
        inh_session{i,1} = 'inh';
    end
end

seshInfo = cell2table([mouse_name, sesh_date, inh_session, sesh_type, correct],...
    'VariableNames', {'mouse_name', 'date', 'inhibition_session', 'session_type', 'correct'});

% Preprocess behavTime
for i = 1: length(behavTime.optoOnset)
   if   isnan(behavTime.optoOnset{i});
       behavTime.optoOnset{i} = 'No_laser';
   elseif behavTime.optoOnset{i} > 0.01;
       behavTime.optoOnset{i} = 0.05;
   else 
       behavTime.optoOnset{i} = 0;
   end
end


% Table to Structure
T = [seshInfo, behavValue, behavTime];
T.Properties.VariableNames = {'mouse_name', 'date', 'inhibition_session', 'session_type', 'correct',...
    'trial_num', 'block_type', 'trial_type', 'response', 'vis_stimulus', 'touch_stimulus',...
    'somOnset', 'somOffset', 'visOnset', 'visOffset', 'stim_onset', 'rLickOnset', 'rLickOffset',...
    'lLickOnset', 'lLickOffset', 'first_lick','opto_stim_onsets'};
S = table2struct(T);

% Save 
dfDir = 'E:\YT080\DF_YT080';
MouseName = se.userData.sessionInfo.MouseName;
seshDate = se.userData.sessionInfo.seshDate;
dfPaths = fullfile(dfDir, ...
    [MouseName,'_', seshDate, '_df.mat']);
save(dfPaths, 'S');

%% Concatenation
MouseName = 'YT080';
MainDir = 'E:\';
dfDir = [MainDir, MouseName, '\DF_', MouseName];
dfFilePaths = MBrowse.Files(dfDir);

df = [];
for i = 1: length(dfFilePaths)
    load(dfFilePaths{i})
    df = [df; S];
    clearvars -except dfFilePaths df MouseName
end

save('E:\YT080\DF_YT080\YT080_df.mat', 'df');