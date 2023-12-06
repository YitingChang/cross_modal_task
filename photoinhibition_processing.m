%%% behavioral performance
% Input: MSessionExplorer of all inhibition sessions for each mouse
% Output: Session-by-session summary of behavioral performance (lick info, trial number, transition info) 
% For each session, behavioral performance must pass the following criteria.
            % Criteria: 
            % Hit rate: lick probability of hit trials
            % Catch rate: lick probability of catch trials (laser only)
            %
            % Without inhibition, 
            % 1. hit rate (visual and tactile stimuli) >= 35% 
            % 2. block performance (touch and light) >= 60% 
            % 3. overall performance >= 65%
            % Without stimuli,
            % 4. catch rate (short and long laser) < hit rate
            

% Directory
Mice = {'YT053' 'YT080' 'YT081' 'YT083' 'YT084' 'YT085' 'YT086'};
Disks = {'D:\' 'E:\' 'E:\' 'F:\' 'F:\' 'E:\' 'E:\'};

for m=1:length(Mice)
    mouseName = Mice{m};
    MainDir = Disks{m};
    seDir = [MainDir, mouseName, '\MSessionExplorer_inhibition'];
%     seFilePaths = MBrowse.Files(seDir);
    cd(seDir)
    seFileLists = struct2cell(dir('*se*'));
    seFileNames = seFileLists(1,:);
    seFilePaths = cellfun(@(x) fullfile(seDir, x),seFileNames,'UniformOutput',false);
    behavSummary = CM.Inhibition.GetBehavInfo(seFilePaths);

    % Save behavior summary
    photoInhPaths = fullfile('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control', ...
        [mouseName,'_photoinhibition']);
    save(photoInhPaths, 'behavSummary');
end

%% Concatenate behavior summary for all mice 
photoInhDir = ['E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control'];
photoInhPaths = MBrowse.Files(photoInhDir);

photoInh_all = table();
for i = 1: length(photoInhPaths)
    load(photoInhPaths{i})
    % Sessions with qualified perfromance are selected
    isPassed = behavSummary.isPassed;
    photoInh_all = [photoInh_all; behavSummary(isPassed,:)];
end
% sortrows(photoInh_all, 'Genotype');
save('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav', 'photoInh_all');