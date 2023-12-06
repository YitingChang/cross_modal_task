%% YT 03/19
    % Simplified versions of  EFcross3_switchArray.m

classdef EFcross3_switchArray_v2 < handle
    
    properties
        mouseName = '';
        sessionDate = '';
        trialNums
        blockType
        trialType
        trialResponse
        visStimType
        somStimType
 
    end
    
    methods (Access = public)
        function obj = EFcross3_switchArray_v2(x, sesDate)
            %
            % function obj = EFcross3Array(x, session_name)
            %
            % Input argument 'x' is either Solo file name string or
            % a structure from loaded Solo file.

            if nargin == 0
                return
            end
            
            if ischar(x)
                x = load(x);
            end
            
            obj.mouseName = x.saved.SavingSection_MouseName;
            obj.sessionDate = sesDate;
            obj.trialNums = x.saved_history.AnalysisSection_NumTrials; %done trials not started trials
            
            blockType = x.saved.TrialTypeSection_previous_block_types; % started trials
            blockType = transpose(blockType);
            blockType(end, :) = []; % changed to done trials: remove the unfinished last trial
            obj.blockType = blockType; 
            
            trialType = x.saved.TrialTypeSection_previous_trial_types; % started trials
            trialType = transpose(trialType); 
            trialType(end,:) = []; % changed to done trials: remove the unfinished last trial
            obj.trialType = trialType;
            
            obj.trialResponse = x.saved.efcross3_switchobj_response_history; % done trials
            obj.visStimType = x.saved_history.TrialTypeSection_VisStimType; % done trials
            obj.somStimType = x.saved_history.TrialTypeSection_StimType; % done trials
            
            %n_trials = length(x.saved_history.ScoringSection_LastTrialEvents);  
        end
    end
end

        