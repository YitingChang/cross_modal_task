classdef Preprocess
    % Modify Many's preprocessing pipeline_YT 022819
    
    methods(Static)
        % Data files to master file
        function intanOps = GetIntanOptions(ksDir)
            
            % customized options for preprocessing Intan Data
            intanOps = MIntan.GetOptions(); % Initialize options
            intanOps(2) = intanOps(1); % Replace aux with another amplifier processing
            intanOps(1).downsampleFactor = 30; % The default Intan sampling rate is 30k Hz and 1000Hz is more than enough for any LFP/EEG analysis 
            intanOps(3).downsampleFactor = 30; % ADC signal: 1000Hz should be enough

            % options for Kilosort
            % 1) Intan amplifier data was acquired with 16-bit resolution but was converted to double precision
            %    numbers. The resulting quantization is 0.195. Kilosort expects 16-bit signed integer. Thus, 
            %    we use the following function to scale amplifier dat
            intanOps(2).isReturn = false; % 2) No need to return a up to integer values and then cast
            %    them to int16. 
            intanOps(2).signalFunc = @(x) int16(x / 0.195); % original amplifier data
            intanOps(2).binFilePath = [ksDir,'\binary_data.dat']; % 3) Specify the path of binary file
        end
        
        function sInfo = GetSessionInfo(MouseName, seshDate, recSite, Genotype, Sex, Manipulation, inhSite, seshType, DoB, behvTask)
            sInfo.MouseName = MouseName;
            sInfo.seshDate = seshDate;
            sInfo.recSite = recSite;
            sInfo.Genotype = Genotype;
            sInfo.Sex = Sex;
            sInfo.Manipulation = Manipulation;
            sInfo.inhSite = inhSite;
            sInfo.seshType = seshType;
            sInfo.DoB = DoB;
            sInfo.behvTask = behvTask;
            
        end
        
        function ADC2SE(intanData, se)

            % Separate channels into different cells
            adcSize = size(intanData.adc_data);
            adcData = mat2cell(intanData.adc_data, adcSize(1), ones(1,adcSize(2)));

            % Derive sampling frequency
            adcFs = 1 / diff(intanData.adc_time(1:2));

            % Find channel names
            adcChanName = {'somStim', 'visStim', 'lLick', 'rLick', 'opto','vTube', 'hTube', 'NaN'}; % Eric's data has 8 channels.
            adcChanName = adcChanName(1:adcSize(2));

            % Import adc time
            adcTime = intanData.adc_time;
            
            % Import trialStartTime
            trialStartTime = se.userData.intanInfo.trialStartTime;
            
            % Make table
            [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(adcTime, adcData, ...
                'DelimiterTimes', trialStartTime,...
                'variableNames', adcChanName);

            % Save to SE
            se.userData.preTaskData.adc = preTb;
            se.SetTable('adc', tb, 'timeSeries');
        end
        
        
        function LFP2SE(intanData, se, sortingMethod)

            % Separate channels into different cells
            ampSize = size(intanData.amplifier_data);
            amplifierData = mat2cell(intanData.amplifier_data, ampSize(1), ones(1,ampSize(2)));

            % Derive sampling frequency
            ampFs = 1 / diff(intanData.amplifier_time(1:2));
            
            switch sortingMethod
                case 'ks'
                    % Find channel names
                    chanMap = se.userData.sessionInfo.channel_map.chanMap;
                    chanMap = num2cell(chanMap);
                    chanMap = cellfun(@(x) strcat('channel_', num2str(x)), chanMap, 'UniformOutput', false);
                    chanMap = chanMap(1:ampSize(2));

                    % Import amplifier time
                    amplifierTime = intanData.amplifier_time;

                    % Import trialStartTime
                    trialStartTime = se.userData.intanInfo.trialStartTime;

                    % Make table
                    [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(...
                        amplifierTime,...
                        amplifierData,...
                        'DelimiterTimes', trialStartTime,...
                        'variableNames', chanMap);
                case 'mclust'
                     % Import amplifier time
                    amplifierTime = intanData.amplifier_time;

                    % Import trialStartTime
                    trialStartTime = se.userData.intanInfo.trialStartTime;

                    % Make table
                    [tb, preTb] = MSessionExplorer.MakeTimeSeriesTable(...
                        amplifierTime,...
                        amplifierData,...
                        'DelimiterTimes', trialStartTime);
            end
            
            % Save to SE
            se.userData.preTaskData.LFP = preTb;
            se.SetTable('LFP', tb, 'timeSeries');
        end
        
        
        
        function Spike2SE(spikeData, se)
            
            % Import trialStartTime
            trialStartTime = se.userData.intanInfo.trialStartTime;
            
            % Creat unit names
            unit_names = num2cell(1:size(spikeData.spike_times, 2));
            unit_names = cellfun(@(x) strcat('unit_', num2str(x)), unit_names, 'UniformOutput', false);
            unit_names = unit_names(1:size(spikeData.spike_times, 2));

            % Make table
            [tb, preTb] = MSessionExplorer.MakeEventTimesTable(...
                spikeData.spike_times, ...
                'delimiterTimes', trialStartTime,...
                'variableNames', unit_names);
            % Save to SE
            se.userData.preTaskData.spikeTime = preTb;
            se.SetTable('spikeTime', tb, 'eventTimes');

        end
        
        
        function BCT2SE(bctData, se)
            % Align bct trial numbers with intan trial numbers
            % Make bev table and save to se
            % bctNums vs. intanNums

            % Bcontrol and Intan data can not be aligned simply by trial
            % number. Need some corrections. (Bcontrol may start from 0 or
            % 1).
            % See TrialArray.m 

            % Import trial numbers from Bcontrol data
            bctNums = bctData.trialNums;
            intNums = transpose(se.userData.intanInfo.trialNums);
            
            % quick fix for sessions with catch trials
            % Repeated trial number for catch trials in Bcontrol
            bctNums = num2cell(0:(length(bctNums)-1)).'; 

            % Convert Trial response from double to cell
            bctTrialResponse = num2cell(bctData.trialResponse);


            % Align bctNums and intNums by trimming front and back trials
            if intNums(1) > bctNums{1}
                frontDiff = abs(intNums(1) - bctNums{1});
                bctNums = bctNums(frontDiff:end);
                bctBlockType = bctData.blockType(frontDiff:end);
                bctTrialType = bctData.trialType(frontDiff:end);
                bctTrialResponse = bctTrialResponse(frontDiff:end);
                bctVisStimType = bctData.visStimType(frontDiff:end);
                bctSomStimType = bctData.somStimType(frontDiff:end);
            elseif intNums(1) < bctNums{1} 
                frontDiff = abs(intNums(1) - bctNums{1});
                bctBlockType = bctData.blockType;
                bctTrialType = bctData.trialType;
                bctTrialResponse = bctTrialResponse;
                bctVisStimType = bctData.visStimType;
                bctSomStimType = bctData.somStimType;
                se.RemoveEpochs(1:frontDiff);
                intNums = intNums(1:frontDiff);
            else
                bctBlockType = bctData.blockType;
                bctTrialType = bctData.trialType;
                bctTrialResponse = bctTrialResponse;
                bctVisStimType = bctData.visStimType;
                bctSomStimType = bctData.somStimType;
            end 

            
            if intNums(end) > bctNums{end}  
                backDiff = abs(intNums(end)-bctNums{end});
                se.RemoveEpochs(length(intNums)-backDiff+2:length(intNums));
                intNums = intNums(1:end-backDiff+1);
            else intNums(end) < bctNums{end}
                backDiff = abs(intNums(end)-bctNums{end});
                bctNums = bctNums(1:end-backDiff-1);
                bctBlockType = bctBlockType(1:end-backDiff-1);
                bctTrialType = bctTrialType(1:end-backDiff-1);
                bctTrialResponse = bctTrialResponse(1:end-backDiff-1);
                bctVisStimType = bctVisStimType(1:end-backDiff-1);
                bctSomStimType = bctSomStimType(1:end-backDiff-1);
            end     


            if size(intNums) == size(bctNums)
                disp('Intan and Bcont trials are correctly aligned')
            else
                error('Intan and Bcont trials are not aligned')
            end

            % Create a table for behavior values
            varNames_bct = {'bct_trialNum', 'blockType', 'trialType', 'response', 'visStimType', 'somStimType'};
            bevTable_bct = table(bctNums, bctBlockType, bctTrialType, bctTrialResponse, bctVisStimType, bctSomStimType,...
                'VariableNames',varNames_bct);

            % update intNums and trialStartTime
            se.userData.intanInfo.trialNums = intNums;
            se.userData.intanInfo.trialStartTime = se.userData.intanInfo.trialStartTime(1:end-backDiff+1);
            
            % Creat a trial type map
            % regular contigency: T-lick Right, V-lick Left
            TTTind = ismember(bevTable_bct.blockType, 'Whisker') & ismember(bevTable_bct.trialType, 'Stim_Som_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [1]) ;

            TTVind = ismember(bevTable_bct.blockType, 'Whisker') & ismember(bevTable_bct.trialType, 'Stim_Som_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [2]) ;

            TTNind = ismember(bevTable_bct.blockType, 'Whisker') & ismember(bevTable_bct.trialType, 'Stim_Som_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [0]) ;

            TVTind = ismember(bevTable_bct.blockType, 'Whisker') & ismember(bevTable_bct.trialType, 'Stim_Vis_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [1]) ;

            TVVind = ismember(bevTable_bct.blockType, 'Whisker') & ismember(bevTable_bct.trialType, 'Stim_Vis_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [2]) ;

            TVNind = ismember(bevTable_bct.blockType, 'Whisker') & ismember(bevTable_bct.trialType, 'Stim_Vis_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [0]) ;

            VTTind = ismember(bevTable_bct.blockType, 'Visual') & ismember(bevTable_bct.trialType, 'Stim_Som_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [1]) ;

            VTVind = ismember(bevTable_bct.blockType, 'Visual') & ismember(bevTable_bct.trialType, 'Stim_Som_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [2]) ;

            VTNind = ismember(bevTable_bct.blockType, 'Visual') & ismember(bevTable_bct.trialType, 'Stim_Som_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [0]) ;

            VVTind = ismember(bevTable_bct.blockType, 'Visual') & ismember(bevTable_bct.trialType, 'Stim_Vis_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [1]) ;

            VVVind = ismember(bevTable_bct.blockType, 'Visual') & ismember(bevTable_bct.trialType, 'Stim_Vis_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [2]) ;

            VVNind = ismember(bevTable_bct.blockType, 'Visual') & ismember(bevTable_bct.trialType, 'Stim_Vis_NoCue') & ...
                ismember(cell2mat(bevTable_bct.response), [0]) ;

            varNames_ttmap = {'TTT','TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                    'VTT','VTV', 'VTN', 'VVT', 'VVV', 'VVN'};

            ttMap = table(TTTind,TTVind, TTNind, TVTind, TVVind, TVNind,...
                    VTTind,VTVind, VTNind, VVTind, VVVind, VVNind, ...
                    'VariableNames',varNames_ttmap); % 12 trial types (2 block types * 2 stim types * 3 responses)
                
            % Merge bct table with ttMap
            bevTable = [bevTable_bct, ttMap];
            
            % Save to SE
            se.SetTable('behavValue', bevTable, 'eventValues');
        end
        
        function GetEventTimes(se)
            % Find event times for visStim, somStim, optoStim, lickOnset
            % Import data
            adcSig = se.GetTable('adc'); % parssed trials
            time = adcSig.time;
            visStim = adcSig.visStim;
            somStim = adcSig.somStim;
            lLick = adcSig.lLick;
            rLick = adcSig.rLick;
            opto = adcSig.opto;

            % Get visual onsets and offsets in each trials
            for k = 1 : length(visStim)
                visInTrial = visStim{k};
                timeInTrial = time{k};
                visOnsetInds{k,1} = find(visInTrial > 0.15, 1);
                % the first one that crosses the threshold of visual stimulus
                visOffsetInds{k,1} = find(visInTrial > 0.15, 1, 'last');
                % the last one that crosses the threshold of visual stimulus

                if isempty(visOnsetInds{k,1}) == 1
                    visOnsetTimes{k,1} = NaN;
                    visOffsetTimes{k,1} = NaN;
                else
                    visOnsetTimes{k,1} = timeInTrial(visOnsetInds{k,1});
                    visOffsetTimes{k,1} = timeInTrial(visOffsetInds{k,1});
                end
            end 
            
            
            % Get som onsets and offsets in each trials
            for k = 1 : length(somStim)
                somInTrial = somStim{k};
                timeInTrial = time{k};
                somOnsetInds{k,1} = find(somInTrial > 0.02, 1);
                somOffsetInds{k,1} = find(somInTrial > 0.02, 1, 'last');

                if isempty(somOnsetInds{k,1}) == 1
                    somOnsetTimes{k,1} = NaN;
                    somOffsetTimes{k,1} = NaN;
                else
                    somOnsetTimes{k,1} = timeInTrial(somOnsetInds{k,1});
                    somOffsetTimes{k,1} = timeInTrial(somOffsetInds{k,1});
                end
            end  
            
            
            % Get stimulus onset
            stimOnsetTimes = {};
            for k = 1 : length(somOnsetTimes)
                if isequaln(somOnsetTimes{k,1},NaN) == 0 && isequaln(visOnsetTimes{k,1}, NaN) == 1  % isequaln treats NaN values as equal,but isequal doesn't
                    stimOnsetTimes{k, 1} = somOnsetTimes{k,1};
                elseif isequaln(somOnsetTimes{k,1}, NaN) == 1 && isequaln(visOnsetTimes{k,1}, NaN) == 0
                    stimOnsetTimes{k, 1} = visOnsetTimes{k,1};
                else
                    stimOnsetTimes{k, 1} = NaN;
                end
            end
            stimOnsetTimes = cell2mat(stimOnsetTimes); % AlignTime function only takes a numeric vactor not a cell array 
            
            % Get right lick onsets and offsets in each trials
            for k = 1 : length(rLick)
                rLickInTrial = rLick{k};
                timeInTrial = time{k};
                rLickThresh = find(rLickInTrial > 0.25);   
                if isempty(rLickThresh) == 1
                    rLickOnsetTimes{k,1} = NaN;
                    rLickOffsetTimes{k,1} = NaN;
                else
                    rLdiff = diff(rLickThresh);
                    rLickOnsetInds = [rLickThresh(1); rLickThresh(find(rLdiff > 80)+1)];
                    rLickOffsetInds = [rLickThresh(find(rLdiff > 80)); rLickThresh(end)];
                    rLickOnsetTimes{k,1} = timeInTrial(rLickOnsetInds);                 
                    rLickOffsetTimes{k,1} = timeInTrial(rLickOffsetInds);
                end
            end  
            
            % Get left lick onsets and offsets in each trials
            for k = 1 : length(lLick)
                lLickInTrial = lLick{k};
                timeInTrial = time{k};
                lLickThresh = find(lLickInTrial > 0.25);
                if isempty(lLickThresh) == 1
                    lLickOnsetTimes{k,1} = NaN;
                    lLickOffsetTimes{k,1} = NaN;
                else
                    lLdiff = diff(lLickThresh);
                    lLickOnsetInds = [lLickThresh(1); lLickThresh(find(lLdiff > 80)+1)];
                    lLickOffsetInds = [lLickThresh(find(lLdiff > 80)); lLickThresh(end)];
                    lLickOnsetTimes{k,1} = timeInTrial(lLickOnsetInds);                 
                    lLickOffsetTimes{k,1} = timeInTrial(lLickOffsetInds);
                end 
            end 
            
            % Get first lick 
            % Select trials based on B-control response 
            % 1) Select correct lick ports
            % 2) Do not select licks after response window 
            %    (0 ~ 2 sec from stimulus onsets)
            % If licks are after response window, response will be 0.

            
            % Load B-control repsonse
            lickResponse = se.GetColumn('behavValue', 'response');
            firstLickTimes = {};
            
            for k = 1 : length(lickResponse)
                if lickResponse{k} == 1 && ~isequaln(rLickOnsetTimes{k}(1), NaN)
                    firstLickInd = find(rLickOnsetTimes{k} > stimOnsetTimes(k), 1); 
                    if ~isempty(firstLickInd) 
                        firstLickTimes{k, 1} = rLickOnsetTimes{k}(firstLickInd);
                    else 
                        firstLickTimes{k, 1} = NaN;
                    end
                elseif lickResponse{k} == 2 && ~isequaln(lLickOnsetTimes{k}(1), NaN) 
                    firstLickInd = find(lLickOnsetTimes{k} > stimOnsetTimes(k), 1); 
                    if ~isempty(firstLickInd) 
                        firstLickTimes{k, 1} = lLickOnsetTimes{k}(firstLickInd);
                    else 
                        firstLickTimes{k, 1} = NaN;
                    end
                else 
                    firstLickTimes{k, 1} = NaN;
                end
            end
            
            firstLickTimes = cell2mat(firstLickTimes); % AlignTime function only takes a numeric vactor not a cell array 
            
            
            % Get opto onsets and offsets 
            for k = 1 : length(opto)
                optoInTrial = opto{k};
                timeInTrial = time{k};
                optoOnsetInds{k,1} = find(optoInTrial > 0.02, 1); % it was 0.1 before 8/9/21
                optoOffsetInds{k,1} = find(optoInTrial > 0.02, 1, 'last');
                
                if isempty(optoOnsetInds{k,1}) == 1
                    optoOnsetTimes{k,1} = NaN;
                    optoOffsetTimes{k,1} = NaN;
                else
                    optoOnsetTimes{k,1} = round(timeInTrial(optoOnsetInds{k,1}),2);
                    optoOffsetTimes{k,1} = round(timeInTrial(optoOffsetInds{k,1}),2);
                end
            end 
            
            % Create behavTime table
            varNames = {'somOnset', 'somOffset', 'visOnset', 'visOffset', 'stimOnset',...
                'rLickOnset', 'rLickOffset', 'lLickOnset', 'lLickOffset', 'firstLick','optoOnset', 'optoOffset'};
            behavTime = table(somOnsetTimes, somOffsetTimes, visOnsetTimes, visOffsetTimes,stimOnsetTimes,...
                rLickOnsetTimes, rLickOffsetTimes, lLickOnsetTimes, lLickOffsetTimes, firstLickTimes, optoOnsetTimes, optoOffsetTimes,...
                'VariableNames',varNames);

            % Add behavTime to se
            se.SetTable('behavTime', behavTime, 'eventTimes');
        end
        
        function TrialQC(se)
            %Quality control of trials
            
            % Remove trials without stimOnset (NaN) due to pre-stim licking
            % within 0.2s before stimulus onset
            stimOnset = se.GetColumn('behavTime', 'stimOnset');
            NaN = isnan(stimOnset);
            se.RemoveEpochs(NaN);
            
            % Remove the last 20 trials (adjust it if drifting happens)
            trialNum = se.GetColumn('behavValue', 'bct_trialNum');
            se.RemoveEpochs(length(trialNum)-20:length(trialNum));
            
        end 
        
    end
end



