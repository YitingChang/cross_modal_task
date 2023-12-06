%% CM_inh_preprocessing
% batch processing

% cd 'F:\Users\Yiting\YitingData\YT086\Intan_YT086'
cd 'F:\YT093\Intan_YT093'

folderNames = dir('*2*');

% % YT053
% inhSite_all = {'sham' 'sham' 'AMM' 'AMM'... % 181115, 16, 18, 19
%             'AMM' 'sham' 'AMM' 'AMM'... % 181120, 21, 22, 23
%             'Prt' 'Prt' 'sham' 'Prt'... % 181124, 26, 27, 28
%             'Prt' 'Prt' 'Prt' 'sham'};  % 181130, 1204, 06, 07

% YT080        
% inhSite_all = {'sham' 'AMM' 'AMM'... % 190503, 05, 06
%             'sham' 'Prt' 'Prt' 'AMM'... % 190507, 08, 09, 10
%             'sham' 'sham' 'AMM' 'Prt'}; % 190511, 12, 13, 14
% 
% % YT081
% inhSite_all = {'sham' 'AMM' 'AMM' 'sham' ... % 190418, 19, 20, 22
%             'Prt' 'Prt' 'sham' 'AMM'... % 190424, 25, 29, 0502
%             'Prt' 'Prt' 'AMM'}; % 190503, 04, 05

% YT084
% % 2019
% inhSite_all = {'sham' 'Prt' 'Prt' 'sham' ... % 191020, 21, 22, 23,
%             'AMM' 'sham' 'MM' 'AMM'... % 191024, 28, 30, 31
%             'Prt' 'sham' 'AMM'... % 191101, 02, 03
%             'MM' 'sham' 'MM' 'S1'... % 191204, 05, 06, 08
%             'S1' 'S1' 'ALM' 'ALM'}; % 191209, 11, 19, 31
% % 2020
% inhSite_all = {'S1' 'ALM' 'ALM'}; % 200103, 06, 08

% YT083 
% 2019
% inhSite_all = {'sham' 'sham' 'S1' 'MM'... % 191205, 06, 07, 08
%             'AMM' 'Prt' 'MM' 'AMM'... % 191209, 10, 11, 12
%             'S1' 'sham' 'Prt' 'Prt'... % 191213, 15, 16, 17
%             'MM' 'AMM' 'sham' 'S1' 'MM'}; % 191218, 19, 22, 26, 31
% 2020
% inhSite_all = {'S1' 'ALM' 'ALM' 'ALM' 'ALM'}; % 200103, 06, 07, 08, 10


% YT085
% inhSite_all = {'sham' 'MM' 'AMM'... %191205, 07, 08
%              'MM' 'sham' 'S1' 'Prt'... % 191210, 11, 12, 13
%              'AMM' 'sham' 'MM'... % 191215, 16, 18
%              'Prt' 'AMM' 'S1' 'S1'... % 191219, 22, 26, 31
%              'MM' 'sham' 'ALM'... % 200102, 03, 06
%              'ALM' 'ALM' 'ALM'}; % 201007, 08, 10

% % YT086
% inhSite_all = {'sham' 'MM' 'AMM'... %191205, 06, 08
%             'MM' 'sham' 'AMM'... %191209, 11, 12
%             'Prt' 'Prt' 'MM'... %191213, 16, 17
%             'sham' 'AMM' 'MM'... % 191218, 19, 31
%             'S1' 'S1' 'ALM'... % 200102, 03, 06
%             'ALM' 'ALM' 'ALM'}; % 201007, 08, 10           


% YT091
% inhSite_all = {'sham' 'sham' 'sham' 'Prt' 'ALM'... % 201124, 25, 26, 27, 28
%     'sham' 'ALM' 'S1' 'Prt'... % 201202, 03, 05, 06
%     'MM' 'AMM' 'S1' 'sham'... % 201207, 08, 09, 11 
%     'ALM' 'AMM' 'MM'}; %201212, 13, 14

% % YT092
% inhSite_all = {'sham' 'sham'... % 201126, 28, 
%      'S1' 'Prt' 'sham' 'S1'... % 201205, 06, 08, 09
%      'ALM' 'AMM' 'MM' 'ALM'... % 201211, 12, 13, 14 
%      'sham' 'wM2'}; % 201217, 19
% 

% YT093
inhSite_all = {'sham' 'sham' 'MM' 'S1'... % 210102, 03, 04, 06
    'AMM' 'sham' 'Prt'}; % 210107, 08, 09
% 200101 (sham) and 200110 (ALM) have a and b sessions

% % YT094
% inhSite_all = {'sham' 'sham' 'sham' 'AMM'... % 201227, 28, 29, 30
%     'MM' 'Prt' 'sham' 'S1'... % 201231, 210102, 03, 06
%     'ALM' 'S1' 'S1' 'Prt'}; % 210108, 09, 13, 14

% % YT095
% inhSite_all = {'sham' 'sham' 'ALM' 'ALM'... %201230, 31, 210101, 02
%     'AMM' 'MM' 'S1' 'MM'... %210103, 04, 06, 07
%     'sham' 'Prt' 'MM' 'Prt'}; %210108, 10, 14, 16 
% % 210109 (AMM) has a and b sessions

%%
% for session = 1:size(folderNames,1)
    % Enter mouse name, session date, recording site, data directory
%     seshDate = folderNames(session).name;
    seshDate = '201230';
    MouseName = 'YT094';
    recSite = 'NaN';
    Genotype = 'PV-Cre-Ai32';
    Sex = 'male';
    Manipulation = 'pre-post';
%     inhSite = inhSite_all{session};
    inhSite = 'AMM';
    seshType = 'inh'; 
    DoB = '200919';
    behvTask = 'tactile detection';
%     behvTask = 'cross-modal selection';

    MainDir = 'F:\';
    IntanDir = [MainDir, MouseName, '\Intan_', MouseName, '\', seshDate];
    BcontDir = [MainDir, MouseName, '\Bcontrol_', MouseName];
    masterDir = [MainDir, MouseName, '\Master_Workspace_', MouseName];
    seDir = [MainDir, MouseName, '\MSessionExplorer_inhibition'];
    ksDir = [MainDir, MouseName, '\Kilosort_', MouseName, '\', seshDate];

    % Import and process Intan data with customized options 
    % customized options for preprocessing Intan Data 
    intanOps = CM.Preprocess.GetIntanOptions(ksDir);
    rhdFileStruct = dir([IntanDir '\*.rhd']);
    rhdFilePaths = arrayfun(@(x) fullfile(x.folder, x.name),rhdFileStruct, 'Uni', false);
    intanStruct = MIntan.ReadRhdFiles(rhdFilePaths, intanOps); 

    %save([masterDir, '\intanStruct.mat'],'intanStruct');

    % Concatenate files to form entire session- masterObj
    % Initialize master file
    masterPath = fullfile(masterDir,[MouseName,'_', seshDate, '_master.mat']);
    masterObj = matfile(masterPath, 'Writable', true);

    % Creat and add sessionInfo
    if isempty(whos(masterObj, 'seshInfo'))
        seshInfo = CM.Preprocess.GetSessionInfo(MouseName, seshDate, recSite, Genotype, Sex, Manipulation, inhSite, seshType, DoB, behvTask);  
        masterObj.seshInfo = seshInfo;
    end

    %Add Intan data
    if isempty(whos(masterObj, 'intan_data')) %&& exist([masterDir, '\intanStruct.mat'], 'file') 
        %intanStruct = fopen([masterDir, '\intanStruct.mat']); 
        masterObj.intan_data = intanStruct;
    end

    % Add Bcontrol data
    if isempty(whos(masterObj, 'bct_data')) 
%         bct_data = CM.EFcross3_switchArray_v2([BcontDir, '\data_@efcross3_switchobj_', MouseName,'_', seshDate, 'a'], seshDate);
        bct_data = CM.EFcross3_switchArray_v2([BcontDir, '\data_@efcross3_switchobj_', MouseName,'_', seshDate], seshDate);
        masterObj.bct_data = bct_data;
        %save([MouseName,'_', sesDate, '_bct.mat'], 'bct_data');
    end


    % Construct SessionExplorers
    % Initialize SE
    cd(seDir);
    se = MSessionExplorer();

    % Load masterObj
    %masterPath = fullfile(masterDir,[MouseName,'_', seshDate, '_master.mat']);
    %masterObj = load(masterPath);

    % Add sessionInfo to se
    se.userData.sessionInfo = masterObj.seshInfo;

    % Load Intan data
    fprintf('Loading intan data\n');
    intanData = masterObj.intan_data;

    % Process delimiter (Dig channel)
    disp('Process delimiter');
    delimiterData = CM.BitCode.ComputeBitCode(intanData.dig_in_data(:,1), 30000); %Should use intanData.info to provide sample rate (30000Hz)
    se.userData.intanInfo.trialStartTime = delimiterData.HexDecOnsetTime;
    se.userData.intanInfo.trialNums = delimiterData.trialNums;
    fprintf('\n');        

    % Process ADC signals
    if ~isempty(intanData.adc_data) && ~ismember('adc', se.tableNames)
        disp('Processing ADC signals');
        CM.Preprocess.ADC2SE(intanData,se);
        fprintf('\n');
    end

    % Process Bcontrol data and create a trial type map
    if ~isempty(masterObj.bct_data) && ~ismember('behavValue', se.tableNames)
        disp('Processing Bcontrol data');
        CM.Preprocess.BCT2SE(masterObj.bct_data, se);
        fprintf('\n');
    end

    % Get event times
    if ~ismember('behavTime', se.tableNames)
        disp('Processing event times');
        CM.Preprocess.GetEventTimes(se);
        fprintf('\n');
    end

    % Set reference time
    se.SetReferenceTime(se.userData.intanInfo.trialStartTime);

    
    % Add stimonset of laser catch trials and catch trials
    % laser catch trials: laser onsets
    % catch trials: trial onsets (0)- tactile detection task only
    % 210820 Stimonsets for catch and laser catch trials were changed due to analyses of pre/post inhibition and censor/grace periods  
    % short laser catch trials: laser offsets
    % long laser catch trials: laser onsets
    % catch trials: trial onsets (0) + 1 - tactile detection task only
    
    stimOnset = se.GetColumn('behavTime', 'stimOnset');
    trialType = se.GetColumn('behavValue', 'trialType');
    optoOnset = se.GetColumn('behavTime', 'optoOnset');
    optoOffset = se.GetColumn('behavTime', 'optoOffset');
    
    new_stimOnset = [];
    for i = 1: length(stimOnset)
        switch behvTask
            case 'tactile detection'
                if strcmp(trialType{i},'catch')|strcmp(trialType{i},'Stim_Vis') % laser catch trial
                    if optoOffset{i} < 1 % short laser catch 
                        new_stimOnset{i,1} = optoOffset{i};
                    else % long laser catch 
                        new_stimOnset{i,1} = optoOnset{i};
                    end
                elseif strcmp(trialType{i},'Stim_Vis_NoCue') % catch trial 
                    new_stimOnset{i,1} = 1;
                else
                    new_stimOnset{i,1} = stimOnset(i);
                end
            case 'cross-modal selection'
                if strcmp(trialType{i},'catch')
                    if optoOffset{i} < 1 % short laser catch 
                        new_stimOnset{i,1} = optoOffset{i};
                    else % long laser catch 
                        new_stimOnset{i,1} = optoOnset{i};
                    end
                else
                    new_stimOnset{i,1} = stimOnset(i);
                end
        end
    end

    se.SetColumn('behavTime', 'stimOnset', cell2mat(new_stimOnset));

    % Quality control of trials (remove trials with pre-stim licking and the
    % last 20 trials)
    CM.Preprocess.TrialQC(se);

    % Align Time to stimulus onset
    se.AlignTime('stimOnset', 'behavTime');


    % Save SE
    sePaths = fullfile(seDir, ...
        [MouseName,'_', seshDate, '_se.mat']);

    save(sePaths, 'se');

    clearvars -except folderNames inhSite_all

% end
