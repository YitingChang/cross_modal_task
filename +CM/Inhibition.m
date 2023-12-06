classdef Inhibition
    methods(Static)
        function behavSummary = GetBehavInfo(seFilePaths)
            % Input: 
            %         seFilePaths: file paths for MSessionExplorers
            % Output: table of behavioral performance 
            %         MouseName: mouse name
            %         seshDate: session date
            %         Genotype: genotype 
            %         inhSite: inhibition site
            %
            %         data: a 4D logical array (b,s,l,o) for lick info in different conditions 
            %               b: respond-to-touch (1) and respond-to-light (2) blocks
            %               s: tactile (1), visual (2), no (3-5) stimuli
            %               l: right (1), left (2), no (3) licks
            %               o: no (1), pre- (2), post- (3) optoinhibition 
            %               For example, data(1,1,1,2) includes all tactile stimulus trials in respond-to-touch blocks and 
            %               optoinhibition is before stimulus onsets. 
            %               One means right lick is detected whereas zero does not for that trial. 
            
            %               Note: For trials without stimuli, 
            %               s=3 for all catch trials, s=4 for short catch trials (0.8 sec), s=5 for long catch trials (2 sec);
            %               o=1 for ITI, o=3 for laser trials  
            %               Response window: 0~2 second from stimulus onsets 
            %
            %         trialNum_tBlock: a 4D array (b,s,l,o) for trial number info in respond-to-touch blocks
            %         trialNum_vBlock: a 4D array (b,s,l,o) for trial number info in respond-to-light blocks
            %         transition_array: a 4D logical array (b,s,l,o) for transition info in different conditions 
            
            % setting
            B=2; % respond-to-touch and respond-to-light blocks
            S=5; % tactile, visual, no stimulus
            L=3; % right, left, no lick
            O=3; % no, pre-, post- inhibition 
            
            for i = 1: length(seFilePaths)
                load(seFilePaths{i})

                %%%%%%%% index for block,stim, opto %%%%%%%%
                isTblock = ismember(se.GetColumn('behavValue', 'blockType'), 'Whisker');
                isVblock = ismember(se.GetColumn('behavValue', 'blockType'), 'Visual');
                block_idx = {isTblock isVblock};
                isTstim = ~isnan(cell2mat(se.GetColumn('behavTime', 'somOnset')));
                isVstim = ~isnan(cell2mat(se.GetColumn('behavTime', 'visOnset'))); 
                isStim = isTstim | isVstim;

                optoOnset = cellfun(@(x) round(x,1), se.GetColumn('behavTime', 'optoOnset'), 'UniformOutput',false); 
                optoOffset = cellfun(@(x) round(x,1), se.GetColumn('behavTime', 'optoOffset'), 'UniformOutput',false); 
                isOpto = ~isnan(cell2mat(optoOnset));
                isPreInh = logical(cell2mat(cellfun(@(x) x<0, optoOnset, 'UniformOutput',false))); 
                isPostInh = logical(cell2mat(cellfun(@(x) x==0, optoOnset, 'UniformOutput',false)));
                opto_idx = {~isOpto isPreInh isPostInh};
                isOptoCatch = isOpto & ~isStim; 
                isShortOptoCatch = isOptoCatch & isPreInh;
                isLongOptoCatch = isOptoCatch & isPostInh;
                stim_idx = {isTstim isVstim isShortOptoCatch isLongOptoCatch};
                
                %%%%%%%% index for lick %%%%%%%%
                % Response from bcontrol data only applied for stimulus trials since different answer windows were used to analyze 
                % laser catch trials (laser but no stimulus). Lick detection from intan data was used for laser catch trials. 
                
                %%% stimulus trials %%%
                response = se.GetColumn('behavValue', 'response');
                isRlick = logical(cell2mat(cellfun(@(x) x==1, response, 'UniformOutput',false))); 
                isLlick = logical(cell2mat(cellfun(@(x) x==2, response, 'UniformOutput',false))); 
                isLick = isRlick | isLlick;
                lick_idx = {isRlick isLlick ~isLick};
                
                %%% laser catch trials %%%
                rLickOnsetTimes = se.GetColumn('behavTime', 'rLickOnset');
                lLickOnsetTimes = se.GetColumn('behavTime', 'lLickOnset');
                isRlick_rw = logical(cell2mat(cellfun(@(x) sum(x<2 & x>0)>0, rLickOnsetTimes, 'UniformOutput',false))); % response window
                isLlick_rw = logical(cell2mat(cellfun(@(x) sum(x<2 & x>0)>0, lLickOnsetTimes, 'UniformOutput',false))); % response window
                isLick_rw = isRlick_rw | isLlick_rw;
                lick_rw_idx = {isRlick_rw isLlick_rw ~isLick_rw};
                
                % short laser trials
                isRlick_iti_short = logical(cell2mat(cellfun(@(x) sum(x<-0.8 & x>-2.8)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isLlick_iti_short = logical(cell2mat(cellfun(@(x) sum(x<-0.8 & x>-2.8)>0, lLickOnsetTimes, 'UniformOutput',false)));
                isLick_iti_short = isRlick_iti_short | isLlick_iti_short;
                lick_iti_short_idx = {isRlick_iti_short isLlick_iti_short ~isLick_iti_short};                
                
                % long laser trials
                isRlick_iti_long = logical(cell2mat(cellfun(@(x) sum(x<0 & x>-2)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isLlick_iti_long = logical(cell2mat(cellfun(@(x) sum(x<0 & x>-2)>0, lLickOnsetTimes, 'UniformOutput',false)));
                isLick_iti_long = isRlick_iti_long | isLlick_iti_long;
                lick_iti_long_idx = {isRlick_iti_long isLlick_iti_long ~isLick_iti_long};
                
                %%%%%%%% lick during censor and grace periods %%%%%%%%
                isRlick_censor = logical(cell2mat(cellfun(@(x) sum(x<0 & x>-0.2)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isLlick_censor = logical(cell2mat(cellfun(@(x) sum(x<0 & x>-0.2)>0, lLickOnsetTimes, 'UniformOutput',false)));
                isRlick_grace = logical(cell2mat(cellfun(@(x) sum(x<0.1 & x>0)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isLlick_grace = logical(cell2mat(cellfun(@(x) sum(x<0.1 & x>0)>0, lLickOnsetTimes, 'UniformOutput',false)));
                isLick_censor = isRlick_censor | isLlick_censor;
                isLick_grace = isRlick_grace | isLlick_grace;
                
                %%%%%%%% Get trial number %%%%%%%%
                % Get trial numbers in each blocks
                BlockChangeInd = ismember(se.GetColumn('behavValue', 'blockType'), 'Change') ;
                tBlockInd = double(isTblock);
                vBlockInd = double(isVblock);

                % Assign block type to the cahnge trials  
                for t = 2: length(tBlockInd)
                    if BlockChangeInd(t) == 1 && tBlockInd(t-1) ==1
                        tBlockInd(t) = 1;
                    elseif BlockChangeInd(t) == 1 && vBlockInd(t-1) ==1
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

                % tHit and vHit
                istHit = isTblock & isTstim & isRlick;
                isvHit = isVblock & isVstim & isLlick;

                % Create trial numbers (after bock switch) and get the first hit trials 
                for tbn = 1:length(tBlock_start_trial)
                    tBlockInd(tBlock_start_trial(tbn):tBlock_end_trial(tbn)) =  cumsum(tBlockInd(tBlock_start_trial(tbn):tBlock_end_trial(tbn)));
                    if ~isempty(find(istHit(tBlock_start_trial(tbn):tBlock_end_trial(tbn)),1))
                        tBlock_start_behv(tbn) = tBlock_start_trial(tbn) + find(istHit(tBlock_start_trial(tbn):tBlock_end_trial(tbn)),1) -1; % first TTT
                    else
                        tBlock_start_behv(tbn) = tBlock_end_trial(tbn)+1; % Mice didn't switch before the end of the block
                    end
                end
                for vbn = 1:length(vBlock_start_trial)
                    vBlockInd(vBlock_start_trial(vbn):vBlock_end_trial(vbn)) =  cumsum(vBlockInd(vBlock_start_trial(vbn):vBlock_end_trial(vbn)));
                    if ~isempty(find(isvHit(vBlock_start_trial(vbn):vBlock_end_trial(vbn)),1))
                        vBlock_start_behv(vbn) = vBlock_start_trial(vbn) + find(isvHit(vBlock_start_trial(vbn):vBlock_end_trial(vbn)),1) -1; % first VVV
                    else
                        vBlock_start_behv(vbn) = vBlock_end_trial(vbn)+1; % Mice didn't switch before the end of the block 
                    end
                end

                % Create transition index for Hit
                transitionInd = repmat(0,length(tBlockInd),1);
                for tbn = 1:length(tBlock_start_behv)
                    if tBlock_start_behv(tbn) ~=1
                        transitionInd(tBlock_start_trial(tbn):tBlock_start_behv(tbn)-1) = 1; 
                        % The first TTT is not included in transition trials.
                    end
                end
                for vbn = 1:length(vBlock_start_behv)
                    if vBlock_start_behv(vbn) ~=1
                        transitionInd(vBlock_start_trial(vbn):vBlock_start_behv(vbn)-1) =  1;
                        % The first VVV is not included in transition trials.
                    end
                end

                %%%%%%%% Get lick result, trial number, transition index %%%%%%%%
                % Trials with stimuli (s=1:2) 
                for b=1:B
                    block = block_idx{b};
                    for s=1:2
                        stim = stim_idx{s};
                        for l=1:L
                            lick = lick_idx{l};
                            for o=1:O
                                opto = opto_idx{o};
                                data(b,s,l,o) = {lick(block & stim & opto & ~isLick_censor & ~isLick_grace)}; 
                                trialNum_tBlock(b,s,l,o) = {tBlockInd(block & stim & lick & opto & ~isLick_censor & ~isLick_grace)};
                                trialNum_vBlock(b,s,l,o) = {vBlockInd(block & stim & lick & opto & ~isLick_censor & ~isLick_grace)};
                                transition_array(b,s,l,o) = {transitionInd(block & stim & lick & opto & ~isLick_censor & ~isLick_grace)};
                            end
                        end
                    end
                end

                % Trials without stimuli (s=3 for short, s=5 for long catch trials); 
                % o=1 for ITI of laser catch trials, 
                % o=2 resonse window after laser offset for short laser catch trials.
                % o=3 resonse window after laser onset for long laser catch trials.
                
                for b=1:B
                    block = block_idx{b};
                    %%% short laser catch 
                    s=3; stim = stim_idx{s};
                    o=1; % iti
                    for l=1:L
                        lick = lick_iti_short_idx{l} ;
                        data(b,s,l,o) = {lick(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_tBlock(b,s,l,o) = {tBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_vBlock(b,s,l,o) = {vBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        transition_array(b,s,l,o) = {transitionInd(block & stim & ~isLick_censor & ~isLick_grace)};
                    end
                    o=2; % response window
                    for l=1:L
                        lick = lick_rw_idx{l} ;
                        data(b,s,l,o) = {lick(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_tBlock(b,s,l,o) = {tBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_vBlock(b,s,l,o) = {vBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        transition_array(b,s,l,o) = {transitionInd(block & stim & ~isLick_censor & ~isLick_grace)};
                    end
                    
                    %%% long laser catch 
                    s=4; stim = stim_idx{s};
                    o=1; % iti
                    for l=1:L
                        lick = lick_iti_long_idx{l} ;
                        data(b,s,l,o) = {lick(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_tBlock(b,s,l,o) = {tBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_vBlock(b,s,l,o) = {vBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        transition_array(b,s,l,o) = {transitionInd(block & stim & ~isLick_censor & ~isLick_grace)};
                    end
                    o=3; % response window
                    for l=1:L
                        lick = lick_rw_idx{l} ;
                        data(b,s,l,o) = {lick(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_tBlock(b,s,l,o) = {tBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        trialNum_vBlock(b,s,l,o) = {vBlockInd(block & stim & ~isLick_censor & ~isLick_grace)};
                        transition_array(b,s,l,o) = {transitionInd(block & stim & ~isLick_censor & ~isLick_grace)};  
                    end
                end
                
                [isAllPassed, behavPerformance] = CM.Inhibition.PerfControl(data);
                behavSummary(i,1) = {se.userData.sessionInfo.MouseName};
                behavSummary(i,2) = {se.userData.sessionInfo.seshDate};
                behavSummary(i,3) = {se.userData.sessionInfo.Genotype};
                behavSummary(i,4) = {se.userData.sessionInfo.inhSite};
                behavSummary(i,5) = {data};
                behavSummary(i,6) = {isAllPassed};
                behavSummary(i,7) = {behavPerformance};
                behavSummary(i,8) = {trialNum_tBlock};
                behavSummary(i,9) = {trialNum_vBlock};
                behavSummary(i,10) = {transition_array};
                clearvars -except seFilePaths behavSummary B S L O
            end

            % Convert to a table
            VarNames = {'MouseName' 'seshDate' 'Genotype' 'inhSite' 'lick' 'isPassed' 'baseline_performance'...
                'trialNum_tBlock' 'trialNum_vBlock' 'transition'};
            behavSummary = cell2table(behavSummary,'VariableNames', VarNames); 
        end
        
        function  [isAllPassed, behavPerformance] = PerfControl(data)
            % Input: a 4D logical array (b,s,l,o) for lick info in different conditions
            % Output: whether the behavioral performance passed the following criteria
            % 
            % Criteria: 
            % Hit rate: lick probability of hit trials
            % Catch rate: lick probability of catch trials (laser only)
            %
            % Without inhibition, 
            % 1. hit rate (visual and tactile stimuli) >= 35% 
            % 2. block performance (touch and light) >=  55% 
            % 3. overall performance >= 60%
            % Without stimuli,
            % 4. catch rate (short and long laser) <= 75% & hit rate
            
            % Threshold for performance
            hitThreshold = 0.35;
            blockThreshold = 0.55;
            overallThreshold = 0.60;
            catchThreshold = 0.75;
            
            % hit rate
            tHit = round(mean(data{1,1,1,1}),2);
            vHit = round(mean(data{2,2,2,1}),2);
            
            % block performance 
            tBlock_perf = round(mean([data{1,1,1,1}; data{1,2,3,1}]),2);
            vBlock_perf = round(mean([data{2,2,2,1}; data{2,1,3,1}]),2);
            
            % overall performance
            overall_perf = round(mean([data{1,1,1,1}; data{1,2,3,1}; data{2,2,2,1}; data{2,1,3,1}]),2); %tHir, vCR, vHit, tCR
            
            % catch rate 
            tCatch_short = round(mean(data{1,3,1,2}),2);
            tCatch_long = round(mean(data{1,4,1,3}),2);
            vCatch_short = round(mean(data{2,3,2,2}),2);
            vCatch_long = round(mean(data{2,4,2,3}),2);
            
            % Performance control
            isPassed(1) = tHit >= hitThreshold;
            isPassed(2) = vHit >= hitThreshold;
            isPassed(3) = tBlock_perf >= blockThreshold;
            isPassed(4) = vBlock_perf >= blockThreshold;
            isPassed(5) = overall_perf >= overallThreshold;
            isPassed(6) = tCatch_short < tHit & tCatch_short <= catchThreshold;
            isPassed(7) = tCatch_long < tHit & tCatch_long <= catchThreshold;
            isPassed(8) = vCatch_short < vHit & vCatch_short <= catchThreshold;
            isPassed(9) = vCatch_long < vHit & vCatch_long <= catchThreshold;
            
            isAllPassed = sum(isPassed) == length(isPassed);
            
            behavPerformance(1) = tHit;
            behavPerformance(2) = vHit;
            behavPerformance(3) = tBlock_perf;
            behavPerformance(4) = vBlock_perf;
            behavPerformance(5) = overall_perf;
            behavPerformance(6) = tCatch_short;
            behavPerformance(7) = tCatch_long;
            behavPerformance(8) = vCatch_short;
            behavPerformance(9) = vCatch_long;
            
        end
        
        
        function behavSummary = GetBehavInfo_tDetection(seFilePaths)
            % Input: 
            %         seFilePaths: file paths for MSessionExplorers
            % Output: 
            %         Table of behavioral performance 
            %         MouseName: mouse name
            %         seshDate: session date
            %         Genotype: genotype 
            %         inhSite: inhibition site
            %
            %         data: a 2D logical array (s,o) for lick info in different conditions 
            %               s: tactile (1), no (2-5) stimuli (2: catch; 3: all laser catch; 4: short laser catch; 5: long laser catch) 
            %               o: no (1), pre- (2), post- (3) optoinhibition 
            %               For example, data(1,2) includes all tactile stimulus trials and 
            %               optoinhibition is before stimulus onsets. 
            %               One means right lick is detected whereas zero does not for that trial. 
            
            %               Notes for lase catch trials:  
            %               s=3 for all laser catch trials, s=4 for short laser catch trials (0.8 sec), 
            %               s=5 for long laser catch trials (2 sec);
            %               o=1 for ITI.  
            %
            %         Stimulus onsets:
            %         stimulus trials: stimulus onsets
            %         catch trials: one second after trial onsets
            %         short laser catch: laser offsets
            %         long laser catch: laser onsets
            %
            %         Response window (relative to stimulus onsets): 
            %         stimulus trials: 0~2 seconds 
            %         catch trials: 0~2 seconds 
            %         laser catch trials: 
            %         o=1: ITI -2~0 seconds from laser onsets
            %         o=3: 0~2 seconds from laser onsets
                        
            % setting
            S=5; 
            O=3; 
            
            for i = 1: length(seFilePaths)
                load(seFilePaths{i})

                %%%%%%%% index for stim,lick,opto %%%%%%%%
                isStim = ~isnan(cell2mat(se.GetColumn('behavTime', 'somOnset')));
                optoOnset = cellfun(@(x) round(x,1), se.GetColumn('behavTime', 'optoOnset'), 'UniformOutput',false); 
                optoOffset = cellfun(@(x) round(x,1), se.GetColumn('behavTime', 'optoOffset'), 'UniformOutput',false); 
                isOpto = ~isnan(cell2mat(optoOnset));
                isPreInh = logical(cell2mat(cellfun(@(x) x<0, optoOnset, 'UniformOutput',false))); 
                isPostInh = logical(cell2mat(cellfun(@(x) x==0, optoOnset, 'UniformOutput',false)));
                isOptoCatch = isOpto & ~isStim;
                isShortOptoCatch = isOptoCatch & isPreInh;
                isLongOptoCatch = isOptoCatch & isPostInh;
                stim_idx = {isStim ~isStim isOptoCatch isShortOptoCatch isLongOptoCatch};
                opto_idx = {~isOpto isPreInh isPostInh}; % for tactile and catch trials
                
                %%%%%%%% lick during different analysis windows %%%%%%%%
                rLickOnsetTimes = se.GetColumn('behavTime', 'rLickOnset');
                isRlick_rw = logical(cell2mat(cellfun(@(x) sum(x<2 & x>0)>0, rLickOnsetTimes, 'UniformOutput',false))); % response window
                isRlick_iti = logical(cell2mat(cellfun(@(x) sum(x<0 & x>-2)>0, rLickOnsetTimes, 'UniformOutput',false))); % intertrial interval
                % short laser trials: ITI and post-laser-onset response window
                isRlick_rw_short = logical(cell2mat(cellfun(@(x) sum(x<1.2 & x>-0.8)>0, rLickOnsetTimes, 'UniformOutput',false))); 
                isRlick_iti_short = logical(cell2mat(cellfun(@(x) sum(x<-0.8 & x>-2.8)>0, rLickOnsetTimes, 'UniformOutput',false)));
                
                %%%%%%%% no lick censor and grace periods %%%%%%%%
                isRlick_censor = logical(cell2mat(cellfun(@(x) sum(x<0 & x>-0.2)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isRlick_grace = logical(cell2mat(cellfun(@(x) sum(x<0.1 & x>0)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isRlick_censor_short = logical(cell2mat(cellfun(@(x) sum(x<-0.8 & x>-1)>0, rLickOnsetTimes, 'UniformOutput',false)));
                isRlick_grace_short = logical(cell2mat(cellfun(@(x) sum(x<-0.7 & x>-0.8)>0, rLickOnsetTimes, 'UniformOutput',false)));  
                
                %%%%%%%% Get lick result %%%%%%%%
                % Tactile and catch trials 
                for s=1:2
                    stim = stim_idx{s};
                    for o=1:O
                        opto = opto_idx{o};
                        data(s,o) = {isRlick_rw(stim & opto & ~isRlick_censor & ~isRlick_grace)}; 
                    end
                end
                
                % Laser catch trials
                % Short laser catch 
                s=4; stim = stim_idx{s};
                data(s,1) = {isRlick_iti_short(stim & ~isRlick_censor_short & ~isRlick_grace_short)};
                data(s,3) = {isRlick_rw_short(stim & ~isRlick_censor_short & ~isRlick_grace_short)};
                % Long laser catch
                s=5; stim = stim_idx{s};
                data(s,1) = {isRlick_iti(stim & ~isRlick_censor & ~isRlick_grace)};
                data(s,3) = {isRlick_rw(stim & ~isRlick_censor & ~isRlick_grace)};
                % All laser catch
                data(3,1) = {[data{4,1}; data{5,1}]};
                data(3,3) = {[data{4,3}; data{5,3}]};
            
                behavPerformance = CM.Inhibition.baseline_performance_tDetection(data);
                isAllPassed = CM.Inhibition.PerfControl_tDetection(data);
                behavSummary(i,1) = {se.userData.sessionInfo.MouseName};
                behavSummary(i,2) = {se.userData.sessionInfo.seshDate};
                behavSummary(i,3) = {se.userData.sessionInfo.Genotype};
                behavSummary(i,4) = {se.userData.sessionInfo.inhSite};
                behavSummary(i,5) = {data};
                behavSummary(i,6) = {isAllPassed};
                behavSummary(i,7) = {behavPerformance};
                
                clearvars -except seFilePaths behavSummary S O
            end

            % Convert to a table
            VarNames = {'MouseName' 'seshDate' 'Genotype' 'inhSite' 'lick_table' 'isPassed' 'baseline_performance'};
            behavSummary = cell2table(behavSummary,'VariableNames', VarNames); 
        end
        
        
        function  isAllPassed = PerfControl_tDetection(data)
            % Input: a 2D logical array (s,o) for lick info in different conditions 
            % Output: whether the behavioral performance passed the following criteria
            % 
            % Criteria: 
            % Hit rate: lick probability of hit trials
            % Catch rate: lick probability of catch trials (no stimulus and laser)
            % Laser catch rate: lick probability of laser catch trials
            
            % Without inhibition, 
            % 1. hit rate >= 35% 
            % 2. overall performance >= 60%
            % Without stimuli,
            % 3. laser catch rate (all, short and long) < hit rate & 75%
            
            % Threshold for performance
            hitThreshold = 0.35;
            overallThreshold = 0.60;
            catchThreshold = 0.75;
            
            % hit rate
            Hit = round(mean(data{1,1}),2);
            
            % false alarm rate
            FA =  round(mean(data{2,1}),2);
            
            % overall performance
            overall_perf = round(mean([data{1,1}; ~data{2,1}]),2); % Hit, CR (~catch)
            
            % laser catch rate 
            optoCatch_all = round(mean(data{3,3}),2); % respond window strating from laser onsets
            optoCatch_short = round(mean(data{4,3}),2); % respond window strating from laser onsets
            optoCatch_long = round(mean(data{5,3}),2); % respond window strating from laser onsets
            optoCatch_short_pre = round(mean(data{2,2}),2); % respond window strating from laser offsets
            
            % Performance control
            isPassed(1) = Hit >= hitThreshold;
            isPassed(2) = overall_perf >= overallThreshold;
            isPassed(3) = optoCatch_all < Hit & optoCatch_all <= catchThreshold;
            isPassed(4) = optoCatch_short < Hit & optoCatch_short <= catchThreshold;
            isPassed(5) = optoCatch_long < Hit & optoCatch_long <= catchThreshold;
            isPassed(6) = optoCatch_short_pre < Hit & optoCatch_short_pre <= catchThreshold;
            
            isAllPassed = sum(isPassed) == length(isPassed);
        end
        
        function behavPerformance = baseline_performance_tDetection(data)

            behavPerformance(1) = round(mean(data{1,1}),2); % hit rate
            behavPerformance(2) = round(mean(~data{2,1}),2); % correct rejection rate (~catch)
            behavPerformance(3) = round(mean([data{1,1}; ~data{2,1}]),2); % overall performance (Hit, CR) 
            behavPerformance(4) = round(mean(data{3,3}),2); % all laser catch rate 
            behavPerformance(5) = round(mean(data{4,3}),2); % short laser catch rate 
            behavPerformance(6) = round(mean(data{5,3}),2); % long laser catch rate 
            behavPerformance(7) = round(mean(data{2,2}),2); % short laser catch rate for pre inhibition

        end
        
    end
end


