%% CM_preprocessing
%% Enter mouse name, session date, recording site, data directory
clear
MouseName = 'EF0148';
seshDate = '190313';
recSite = 'left POm';
Genotype = 'Emx-Cre;Ai32';
Sex = 'male';
Manipulation = 'NaN';
inhSite= 'NaN';

seshType= 'NaN';
DoB = '180627';


MainDir = 'E:\';
IntanDir = [MainDir, MouseName, '\Intan_', MouseName, '\', seshDate];
ksDir = [MainDir, MouseName, '\Kilosort_', MouseName, '\', seshDate];
BcontDir = [MainDir, MouseName, '\Bcontrol_', MouseName];
masterDir = [MainDir, MouseName, '\Master_Workspace_', MouseName];
seDir = [MainDir, MouseName, '\MSessionExplorer_', MouseName];

% Make subfolders
mkdir(ksDir);
mkdir(masterDir);

%% Import and process Intan data with customized options 
intanOps = CM.Preprocess.GetIntanOptions(ksDir);
rhdFilePaths = MBrowse.Files(IntanDir);
intanStruct = MIntan.ReadRhdFiles(rhdFilePaths, intanOps); 

%save([masterDir, '\intanStruct.mat'],'intanStruct'); 
%% Kilosort
 
MKilosort.Sort([ksDir,'\binary_data.dat']); % select binary_data.dat  
% Kilosort outputs:
% chan map,channel positions, spike templates, spike times,amplitude etc.  
 
% Manual processing for Kilosort   
% Run this part in Phy
% Open terminal and activate phy environment  
% activate phy
% change the drive in command prompt  F:
% cd F:\YT086\Kilosort_YT086\200318 (KsDir)
% phy template-gui params.py

%% Concatenate files to form entire session- masterObj

% Initialize master file
masterPath = fullfile(masterDir,[MouseName,'_', seshDate, '_master.mat']);
masterObj = matfile(masterPath, 'Writable', true);

% Creat session Information 
seshInfo = CM.Preprocess.GetSessionInfo(MouseName, seshDate, recSite, Genotype, Sex, Manipulation, inhSite, seshType, DoB);  
seshInfo.channel_map = load([ksDir,'\chanMap.mat']); 

%Add Intan data
if isempty(whos(masterObj, 'intan_data')) %&& exist([masterDir, '\intanStruct.mat'], 'file') 
    %intanStruct = load([masterDir, '\intanStruct.mat']); 
    masterObj.intan_data = intanStruct;
end

% Add Kilosort and TemplateGUI outputs
if isempty(whos(masterObj, 'spike_data')) 
    spikeData = MKilosort.ImportResults(ksDir); % Import the manual sorting result
    masterObj.spike_data = spikeData;
end
 
% Add session Information
if isempty(whos(masterObj, 'seshInfo'))
    seshInfo.unit_channel_id = spikeData.info.unit_channel_ind; %Indices of primary channel (after mapping) for each unit.
    seshInfo.unit_position = (seshInfo.unit_channel_id - 1)*20; % unit positions  
    masterObj.seshInfo = seshInfo;
end

% Add Bcontrol data
if isempty(whos(masterObj, 'bct_data')) 
    bct_data = CM.EFcross3_switchArray_v2([BcontDir, '\data_@efcross3_switchobj_', MouseName,'_', seshDate, 'a'], seshDate);
    masterObj.bct_data = bct_data;
end
%% Construct SessionExplorers

% Initialize SE
cd(seDir);
se = MSessionExplorer();

% Load masterObj
% masterObj = load(masterPath);

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


% Process LFP signals
if ~isempty(intanData.amplifier_data) && ~ismember('LFP', se.tableNames)
    disp('Processing LFP signals');
    CM.Preprocess.LFP2SE(intanData, se, 'ks');
    fprintf('\n');
end

   
% Load and process spike data
if ~isempty(masterObj.spike_data) && ~ismember('spikeTime', se.tableNames)
    disp('Loading spike data');
    spikeData = masterObj.spike_data;
    disp('Processing spike data');
    CM.Preprocess.Spike2SE(spikeData, se);
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

% Quality control of trials (remove trials with pre-stim licking and the
% last 20 trials)
CM.Preprocess.TrialQC(se);

% Align Time to stimulus onset
se.AlignTime('stimOnset', 'behavTime');

%Remove trials if drifting happens
%  r = ; % the last r trials to be removed 
%  trialNum = se.GetColumn('behavValue', 'bct_trialNum');
%  se.RemoveEpochs(length(trialNum)-r:length(trialNum));
         
%% Save SE
sePaths = fullfile(seDir, ...
    [MouseName,'_', seshDate, '_se.mat']); 
 
save(sePaths, 'se');
