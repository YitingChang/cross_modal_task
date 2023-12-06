%% Silicon probe recordings
% File paths for MSessionExplorers
load('E:\seFilePaths_SiProbe'); 

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% load('E:\Data_array\AllTrials_100bin_noSmooth\data_array.mat')

% Select sessions with interested recording areas
recSites = {'left S1', 'left S2', 'left wM2', 'left ALM'};
u = data_all(ismember(data_all.recSite, recSites),:);

% Get information for each mouse
MouseNames = unique(u.MouseName);
for m = 1:length(MouseNames)
    MouseName = MouseNames{m,1};
    % recording sessions
    v = u(strcmp(u.MouseName, MouseName),:);
    seshDates = v.seshDate;
    if size(seshDates,1) == 1
        seshDate_range = seshDates;
    else
        seshDate_first = seshDates{1,1};
        seshDate_last = seshDates{end,1};
        seshDate_range = [seshDate_first ' - ' seshDate_last];
    end
    
    % DoB, sex, genotype
    seFilePaths = seFilePaths_all(strcmp(seFilePaths_all.MouseName, MouseName),3); 
    seFilePath = seFilePaths{1,1}{1,1}{1,1};
    load(seFilePath)
    genotype = se.userData.sessionInfo.Genotype;
    DoB = se.userData.sessionInfo.DoB;
    sex = se.userData.sessionInfo.sex;
    summary{m,1} = MouseName;
    summary{m,2} = sex;
    summary{m,3} = genotype;
    summary{m,4} = DoB;
    summary{m,5} = seshDate_range;
    summary{m,6} = seshDates;
    
end

VarNames = {'MouseName' 'Sex' 'Genotype' 'DoB' 'seshDate_range' 'seshDates'};
summary = cell2table(summary,'VariableNames', VarNames); 
save('E:\ephy_mouseInfo', 'summary')

mouseInfo = summary(:,1:5);
writetable(mouseInfo,'E:\ephy_mouseInfo.xlsx')

%% Iinhibition - cross-modal selection 

load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav')
MouseNames = unique(photoInh_all.MouseName);
for m = 1:length(MouseNames)
    MouseName = MouseNames{m,1};
    v = photoInh_all(strcmp(photoInh_all.MouseName, MouseName),:);
    genotype = v.Genotype{1,1};
    seshDates = v.seshDate;
    if size(seshDates,1) == 1
        seshDate_range = seshDates;
    else
        seshDate_first = seshDates{1,1};
        seshDate_last = seshDates{end,1};
        seshDate_range = [seshDate_first ' - ' seshDate_last];
    end
    
    % DoB, sex, genotype
    summary{m,1} = MouseName;
    summary{m,3} = genotype;
    summary{m,5} = seshDate_range;
    summary{m,6} = seshDates;
end

VarNames = {'MouseName' 'Sex' 'Genotype' 'DoB' 'seshDate_range' 'seshDates'};
summary = cell2table(summary,'VariableNames', VarNames); 
save('E:\CM_inh_mouseInfo', 'summary')

mouseInfo = summary(:,1:5);
writetable(mouseInfo,'E:\CM_inh_mouseInfo.xlsx')
%% Inhibition - tactile detection 
load('E:\CM_Photoinhibition_Analysis\Tactile_detection\photoInh_behav')
MouseNames = unique(photoInh_all.MouseName);
for m = 1:length(MouseNames)
    MouseName = MouseNames{m,1};
    v = photoInh_all(strcmp(photoInh_all.MouseName, MouseName),:);
    genotype = v.Genotype{1,1};
    seshDates = v.seshDate;
    if size(seshDates,1) == 1
        seshDate_range = seshDates;
    else
        seshDate_first = seshDates{1,1};
        seshDate_last = seshDates{end,1};
        seshDate_range = [seshDate_first ' - ' seshDate_last];
    end
    
    % DoB, sex, genotype
    summary{m,1} = MouseName;
    summary{m,3} = genotype;
    summary{m,5} = seshDate_range;
    summary{m,6} = seshDates;
end

VarNames = {'MouseName' 'Sex' 'Genotype' 'DoB' 'seshDate_range' 'seshDates'};
summary = cell2table(summary,'VariableNames', VarNames); 
save('E:\Tac_inh_mouseInfo', 'summary')

mouseInfo = summary(:,1:5);
writetable(mouseInfo,'E:\Tac_inh_mouseInfo.xlsx')


        