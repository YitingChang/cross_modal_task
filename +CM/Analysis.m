%% Create a class
classdef Analysis
    
    methods(Static)
        function spkRs = GetSpikeRate(se, bin, tStart, tEnd)
            
            % Inputs
            % bin size (sec), time window of interest (from tStart to tEnd)
            %
            % Outputs
            % spike rate (spkRs) 
            
            % Resample
            tEdges = tStart: bin : tEnd;
            spkRs = se.ResampleEventTimes('spikeTime', tEdges);
        end

        
        function normalized_spkRs = zScore(se, bin, spkRs)
            % Normalization
            % Calculate baseline for each units (average all trial types-> The weights
            % of differnt trial types differ. larger in hit and cr trials)
            % Output
            % Normalized spike rate (zScore)
            
            tEdges_baseline = -1 : bin : 0; % baseline -1~ 0 sec % it was -.5~ 0 sec before 9/26/19
            spkRs_baseline = se.ResampleEventTimes('spikeTime', tEdges_baseline);

            bl = [];
            for k = 2 : width(spkRs_baseline)
                baseline = cellfun(@mean,spkRs_baseline{:,k});
                mean_baseline = mean(baseline);
                std_baseline = std(baseline);
                bl{1,k} = mean_baseline;
                bl{2,k} = std_baseline;
            end

            normalized_spkRs = [spkRs.time];
            for k = 2 : width(spkRs)
                normalized_unit = cellfun(@(x) (x -bl{1,k})/ bl{2, k}, spkRs{:,k}, 'UniformOutput', false);
                normalized_spkRs = [normalized_spkRs normalized_unit];
            end
            normalized_spkRs = cell2table(normalized_spkRs,'VariableNames', spkRs.Properties.VariableNames); 
            normalized_spkRs = normalized_spkRs(:,2:end); % remove time
        end
        
        function spkRate_tt = GetUnitSummary(se, spkRate, varargin) 
            % Input
            % spike data (spkData) can be spike rates or normalized spike
            % rates (zScore)
            % varargin: use bootsrapping or not (Improvement: shorten time
            % window and the number of trial types, choose units, and use 
            % parallel processing to reduce processing time)
       
            % Output
            % a table includes mean, std, sem, the lower bound of sem (lbSem)
            % and the higher bound of sem (hbSem)of all trial types for
            % each unit
            
            % Load trial type index
            trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                  'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
            ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
                'VariableNames', trialTypes);
            
            for tt = 1:size(ttInd, 2)
                ind = ttInd{:, tt};
                mean_tt = [];
                std_tt = [];
                lbCI_tt = [];
                hbCI_tt = [];
                sig_ind_tt = [];

                for unit = 1 : size(spkRate, 2) % The first column is time.
                    s = spkRate{:,unit}; % Select unit
                    s = s(ind); % Select trial types

                    unit_tb = [];

                    for trial = 1 : length(s)
                        unit_tb = [unit_tb, s{trial,1}];
                    end
                    
                    switch nargin
                        % Get mean, std, sem of all trial types for each unit
                        case 2                      
                            mean_tt_unit = mean(unit_tb, 2); % mean trial-type trace for each unit
                            std_tt_unit = std(unit_tb,0,2); % std trial-type trace for each unit
                            CI_tt_unit = 1.96*std_tt_unit/ sqrt(length(unit_tb)); % changed from sem to 95% CI
                            %sem_tt_unit = std_tt_unit/ sqrt(length(unit_tb));
                            lbCI_tt_unit = mean_tt_unit - CI_tt_unit;
                            hbCI_tt_unit = mean_tt_unit + CI_tt_unit;
                            sig_ind_tt_unit = lbCI_tt_unit.*hbCI_tt_unit >0; % if 95%CI contains 0, lbCI*hbCi will be negative
                            
                        % Get mean and 95% CI by bootstrapping for each time bin  
                        case 3
                            rng(22); % control random number generation
                            nboot = 1000; % nboot
                            unit_tb2 = num2cell(unit_tb,2);
                            bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), unit_tb2, 'UniformOutput', false);
                            bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
                            mean_tt_unit =  cellfun(@(x) mean(x), bootstat_sorted,'UniformOutput', false);
                            std_tt_unit = cellfun(@(x) std(x), bootstat_sorted,'UniformOutput', false);
                            lbCI_tt_unit = cellfun(@(x)x(nboot*0.025), bootstat_sorted,'UniformOutput', false);
                            hbCI_tt_unit = cellfun(@(x)x(nboot*0.975), bootstat_sorted,'UniformOutput', false);
                            sig_ind_tt_unit = cell2mat(lbCI_tt_unit).*cell2mat(hbCI_tt_unit) >0; % if 95%CI contains 0, lbCI*hbCi will be negative                          
                    end
                    % Add each unit trace to a table 
                    mean_tt = [mean_tt, mean_tt_unit];
                    std_tt = [std_tt, std_tt_unit];
                    lbCI_tt = [lbCI_tt, lbCI_tt_unit];
                    hbCI_tt = [hbCI_tt, hbCI_tt_unit];
                    sig_ind_tt = [sig_ind_tt, sig_ind_tt_unit];                  
                    
                end

                spkRate_tt{1, tt} = mean_tt;
                spkRate_tt{2, tt} = std_tt;
                spkRate_tt{3, tt} = lbCI_tt;
                spkRate_tt{4, tt} = hbCI_tt;
                spkRate_tt{5, tt} = sig_ind_tt;
            end

            % Create a table
            spkRate_tt = cell2table(spkRate_tt, 'VariableNames', trialTypes, 'RowNames', {'mean','std', 'lb95CI', 'hb95CI', 'sig_ind'}); 
        end
        
        function isResponsive = GetResponsiveUnit(se, timeWindow) 
            % Determine responsive units
            % more than 5 trials with maximal firing rate > 5 Hz
            spkRs_isResponsive = CM.Analysis.GetSpikeRate(se, diff(timeWindow), timeWindow(1), timeWindow(2)); 
            for k = 1 : width(spkRs_isResponsive)-1
                a = find(cell2mat(spkRs_isResponsive{:,k+1}) > 5*diff(timeWindow));
                isResponsive(k) = length(a)>5;
            end 
        end
        
        
        function zScore_indConcat(MouseName, seDir, bin, timeWindow, varargin)

            % Inputs
            %   MouseName    Mouse name
            %   seDir        Directory of session explore
            %   bin          Bin size (sec)
            %   timeWindow   The period of time for analysis (sec) [tStart, tEnd]
            %   varargin     Align to lick onset or not
            % Outputs
            %   Concatenation of zScores of all responsive units from individual mice
            %   zScore results include mean, std, sem, the lower bound of sem (lbSem)
            %   and the higher bound of sem (hbSem)of all trial types

            % Concatenation
            % Load seFilePaths
            seFilePaths = MBrowse.Files(seDir);

            zScore = [];
            for i = 1: length(seFilePaths)
                load(seFilePaths{i})

                switch nargin
                    case 4

                        spkRs = CM.Analysis.GetSpikeRate(se, bin, timeWindow(1), timeWindow(2)); % It was [-.5 2.5] before 9/26/19. 
                        normalized_spkRs = CM.Analysis.zScore(se, bin, spkRs);

                    case 5 % Align to lick onsets
                        % Get baseline (0.5s before stim onset) 
                        tEdges_baseline = -1 : bin : 0; % baseline -1~ 0 sec
                        spkRs_baseline = se.ResampleEventTimes('spikeTime', tEdges_baseline);

                        bl = [];
                        for k = 2 : width(spkRs_baseline)
                            baseline = cellfun(@mean,spkRs_baseline{:,k});
                            mean_baseline = mean(baseline);
                            std_baseline = std(baseline);
                            bl{1,k} = mean_baseline;
                            bl{2,k} = std_baseline;
                        end
                        %%%%%%%%% Align to lick onsets %%%%%%%%%%%%
                        firstLick = se.GetColumn('behavTime', 'firstLick'); % replace NaN with 0 for no lick trials
                        firstLick(isnan(firstLick)) = 0;
                        se.AlignTime(firstLick);

                        % Get spike rate
                        spkRs = CM.Analysis.GetSpikeRate(se, bin,  timeWindow(1), timeWindow(2));

                        % Normalization
                        normalized_spkRs = [spkRs.time];
                        for k = 2 : width(spkRs)
                            normalized_unit = cellfun(@(x) (x -bl{1,k})/ bl{2, k}, spkRs{:,k}, 'UniformOutput', false);
                            normalized_spkRs = [normalized_spkRs normalized_unit];
                        end
                        normalized_spkRs = cell2table(normalized_spkRs,'VariableNames', spkRs.Properties.VariableNames); 
                        normalized_spkRs = normalized_spkRs(:,2:end); % remove time                     
                end


                % Get zScore summary
                spkRate_tt = CM.Analysis.GetUnitSummary(se, normalized_spkRs);
                
                % Determine whether it is responsive unit (>= 5Hz at least
                % for five trials)
                isResponsive = CM.Analysis.GetResponsiveUnit(se, timeWindow); 
                
                % Load session information 
                MouseName = se.userData.sessionInfo.MouseName;
                seshDate = se.userData.sessionInfo.seshDate;
                recSite = se.userData.sessionInfo.recSite;

                % Combine with zScore
                new_zScore = [{MouseName} {seshDate} {recSite} table2cell(spkRate_tt(1,:)) {isResponsive}];

                % Concatenation 
                zScore = [zScore; new_zScore];

                clearvars -except seFilePaths zScore MouseName bin timeWindow          
            end
            % Convert to a table
            trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                              'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
            VarNames = ['MouseName' 'seshDate' 'recSite' trialTypes 'isResponsive'];
            zScore = cell2table(zScore,'VariableNames', VarNames);

            % Save zScore summary for each mouse 
            switch nargin
                case 4
                    zScorePaths = fullfile('E:\CM_NeuralActivity_Analysis\zScore_responsive', ...
                        [MouseName,'_zScore']);
                    save(zScorePaths, 'zScore');
                case 5
                    zScorePaths = fullfile('E:\CM_NeuralActivity_Analysis\zScore_responsive', ...
                        [MouseName,'_zScore_lickOnset']);
                    save(zScorePaths, 'zScore');
            end
                    
        end
  
            
        function zScore_allConcat (zScoreDir, varargin)
            zScorePaths = MBrowse.Files(zScoreDir);
            zScore_all = table();
            for i = 1: length(zScorePaths)
                load(zScorePaths{i})
                zScore_all = [zScore_all; zScore];
            end
            % Save zScore summary for all mice 
            switch nargin
                case 1
                    zScoreSavePath = fullfile(zScoreDir, 'zScore_SiProbe');
                    save(zScoreSavePath, 'zScore_all');
                case 2
                    zScoreSavePath = fullfile(zScoreDir, 'zScore_lickOnset');
                    save(zScoreSavePath, 'zScore_all');
            end
        end
        
        function ranSelect(zScore_all, n)

            % Randomly select n number of units from tetrode dataset 
            % Load trial type index
            trialTypes = zScore_all.Properties.VariableNames(4:end);
            
            for tt = 1:size(trialTypes, 2)               
                % Concatenation
                zScore_tetrode_tt = [];
                for j = 1:size(zScore_all,1) % loop each session
                    % get the number of units in each session
                    n_units = size(zScore_all{j,4}{1,1},2);
                    if ~isempty(zScore_all{j,tt+3}{1,1})
                        zScore_tetrode_tt = [zScore_tetrode_tt zScore_all{j,tt+3}{1,1}];
                    else
                        zScore_tetrode_tt = [zScore_tetrode_tt NaN(350,n_units)]; % assign NaN to sessions without certain trial types
                    end
                end
                zScore_tetrode{tt} = zScore_tetrode_tt;
            end
            
            rng(19); % control random number generation
            randomInd = randperm(size(zScore_tetrode{1,1},2),n); % randomly select n number from 1 tototal number of units
            for tt = 1:size(trialTypes, 2)  
                if tt == 13
                    zScore_tetrode_selected{tt} = logical(zScore_tetrode{tt}(:,randomInd));
                else
                    zScore_tetrode_selected{tt} = zScore_tetrode{tt}(:,randomInd);
                end
            end
            recSite = zScore_all.recSite{1};
            zScore_tetrode_selected_tb = [{'Random'} {'Random'} {recSite} zScore_tetrode_selected]; 
            % Convert to a table
            VarNames = {'MouseName', 'seshDate', 'recSite', trialTypes{1,:}};
            zScore = cell2table(zScore_tetrode_selected_tb,'VariableNames', VarNames);

            % Save zScore summary for tetrode recordings
            zScorePaths = fullfile('E:\CM_NeuralActivity_Analysis\zScore_responsive',...
                ['Random', '_', recSite, '_zScore']);
            save(zScorePaths, 'zScore');              
        end
        
        function tetrode2siConcat(zScoreDir, zScore_SiProbe)
            % Input
            %   zScoreDir: where to find block type index of tetrode dataset
            %   zScore_SiProbe: block type index of silicone probe dataset
            % Concatenation 
            zScorePaths = MBrowse.Files(zScoreDir);
            zScore_all = zScore_SiProbe;
            for i = 1: length(zScorePaths)
                load(zScorePaths{i})
                zScore_all = [zScore_all; zScore];
            end

            % Save btInx summary for all mice 
            save('E:\CM_NeuralActivity_Analysis\zScore_responsive\zScore_mixed', 'zScore_all');
        end
        
        function zScore_recSite = GetGrandAvg(zScore_all, recSite)
            warning('off')
            for i = 1: length(recSite)
                u = zScore_all(strcmp(zScore_all.recSite, recSite{i}), :);
                % Get mean, std, sem for all trial types
                % Load trial types
                trialTypes = zScore_all.Properties.VariableNames(4:15);

                for tt = 1: length(trialTypes) 
                    % Concatenation 
                    v = [];
                    for session = 1:height(u)
                        % Select responsive units
                        isResponsive = u{session, 16}{1,1};
                        new_v = u{session, 3 + tt}{1,1}(:,isResponsive); % TrialTypes start from the 4th column because the first three columns are MouseName, seshDate, recSite
                        if sum(isnan(new_v),'all') ~= 0
                        % remove units that do not have tt trial types
                        % (NaN)
                            n_timebin = size(new_v,1);
                            new_v = new_v(~isnan(new_v));
                            new_v = reshape(new_v, [n_timebin, length(new_v)/n_timebin]);
                        end
                        v = [v new_v];
                    end
                    
                    % bootstrapping for each time bin 
                    rng(17); % control random number generation
                    nboot = 1000; % nboot
                    v2 = num2cell(v,2);
                    bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), v2, 'UniformOutput', false);
                    bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
                     bootstrap_summary = cellfun(@(x) [mean(x), x(nboot*0.025), x(nboot*0.975)], bootstat_sorted,'UniformOutput', false);
                    bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI           
            
                   
                    % Get zScore summary
                    zScore_tt{1, tt} = bootstrap_summary(:,1); % mean
                    zScore_tt{2, tt} = bootstrap_summary(:,2); % lower CI
                    zScore_tt{3, tt} = bootstrap_summary(:,3);; % upper CI             

%                     zScore_tt{1, tt} = mean(v, 2); % mean
%                     zScore_tt{2, tt} = std(v,0,2); % std
%                     zScore_tt{3, tt} = std(v,0,2)/ sqrt(length(v)); %sem
%                     zScore_tt{4, tt} = mean(v, 2) - std(v,0,2)/ sqrt(length(v)); %lower bound of sem
%                     zScore_tt{5, tt} = mean(v, 2) + std(v,0,2)/ sqrt(length(v)); % higher bound of sem
                end
                % Create a table
                zScore_tt = cell2table(zScore_tt, 'VariableNames', trialTypes,...
                    'RowNames', {'mean','lower CI', 'upper CI'});

                % Add zScore summary of each recording sites
                zScore_recSite{i,1} = recSite{i};
                zScore_recSite{i,2} = zScore_tt;    
                clearvars -except zScore_all recSite zScore_recSite 
            end
            warning('on') 
            
        end
            

        function plot_correctTrials (zScore_recSite, type, timeWindow, save_path)
                % Inputs
                %   zScore summary for each recording site
                %   plot based on stimulus types or block types: stimulus or block  
                %   time window (sec) ex: [-0.25 1] 0: stimulus onset
                %   save_path ex:
                %   'E:\CM_NeuralActivity_Analysis\zScore\GrandAvg'
                % Outputs
                %   figure saved in a specific path
            
            time = -1:0.01:2.49;
            for i = 1: size(zScore_recSite,1)
                zScore_tt = zScore_recSite{i,2};
                recSite = zScore_recSite{i,1};
                figure('Position', [0,0, 800, 500]); % set figure size
                switch type
                    case 'block'
                        subplot(1,2,1)
                        plot(time, zScore_tt.TTT{1}, 'b',time, zScore_tt.TVN{1},'r:');
                        hold on
                        sem_TTT = [zScore_tt.TTT{2}.', fliplr(zScore_tt.TTT{3}.')];
                        sem_TVN = [zScore_tt.TVN{2}.', fliplr(zScore_tt.TVN{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_TTT, 'b','FaceAlpha', 0.3, 'Linestyle', 'none');
                        fill(time_bw, sem_TVN, 'r','FaceAlpha', 0.1, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Touch Block', 'FontSize',16);
                        hold off

                        subplot(1,2,2)
                        plot(time, zScore_tt.VVV{1},'r',time,zScore_tt.VTN{1},'b:');
                        hold on
                        sem_VVV = [zScore_tt.VVV{2}.', fliplr(zScore_tt.VVV{3}.')];
                        sem_VTN = [zScore_tt.VTN{2}.', fliplr(zScore_tt.VTN{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_VVV, 'r', 'FaceAlpha', 0.3, 'Linestyle', 'none');
                        fill(time_bw, sem_VTN, 'b','FaceAlpha', 0.1, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Light Block', 'FontSize',16);
                        hold off

                        sgtitle(recSite, 'FontSize',20);

                        % Save fig
                        zScoreFigPath = fullfile(save_path,[recSite,'_CorrectTrials_BlockType.pdf']);
                        print(zScoreFigPath,'-dpdf','-painters','-loose');
                    case 'stimulus'
                        subplot(1,2,1)
                        plot(time, zScore_tt.TTT{1}, 'b',time,zScore_tt.VTN{1},'b:');
                        hold on
                        sem_TTT = [zScore_tt.TTT{2}.', fliplr(zScore_tt.TTT{3}.')];
                        sem_VTN = [zScore_tt.VTN{2}.', fliplr(zScore_tt.VTN{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_TTT, 'b','FaceAlpha', 0.3, 'Linestyle', 'none');
                        fill(time_bw, sem_VTN, 'b','FaceAlpha', 0.1, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Tactile Stim', 'FontSize',16);
                        hold off

                        subplot(1,2,2)
                        plot(time, zScore_tt.VVV{1},'r',time, zScore_tt.TVN{1},'r:');
                        hold on
                        sem_VVV = [zScore_tt.VVV{2}.', fliplr(zScore_tt.VVV{3}.')];
                        sem_TVN = [zScore_tt.TVN{2}.', fliplr(zScore_tt.TVN{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_VVV, 'r','FaceAlpha', 0.3, 'Linestyle', 'none');
                        fill(time_bw, sem_TVN, 'r','FaceAlpha', 0.1, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Visual Stim', 'FontSize',16);
                        hold off

                        sgtitle(recSite, 'FontSize',20);

                        % Save fig
                        zScoreFigPath = fullfile(save_path,[recSite,'_CorrectTrials_StimType.pdf']);
                        print(zScoreFigPath,'-dpdf','-painters','-loose');
                end
            end    
        end
        
        function plot_allTrials (zScore_recSite, type, timeWindow, save_path)
            time = -1:0.01:2.49;
            for i = 1: size(zScore_recSite,1)
                zScore_tt = zScore_recSite{i,2};
                recSite = zScore_recSite{i,1};
                figure('Position', [0,0, 800, 500]); % set figure size
                switch type
                    case 'block'
                        subplot(1,2,1)
                        plot(time, zScore_tt.TTT{1}, 'b', time, zScore_tt.TTN{1},'k',...
                            time,zScore_tt.TVN{1},'r', time, zScore_tt.TVT{1},'g');
                        hold on
                        sem_TTT = [zScore_tt.TTT{2}.', fliplr(zScore_tt.TTT{3}.')];
                        sem_TTN = [zScore_tt.TTN{2}.', fliplr(zScore_tt.TTN{3}.')];
                        sem_TVN = [zScore_tt.TVN{2}.', fliplr(zScore_tt.TVN{3}.')];
                        sem_TVT = [zScore_tt.TVT{2}.', fliplr(zScore_tt.TVT{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_TTT, 'b', time_bw, sem_TTN, 'k',...
                            time_bw, sem_TVN, 'r',time_bw, sem_TVT, 'g',...
                            'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Touch Block', 'FontSize',16);
                        hold off

                        subplot(1,2,2)
                        plot(time, zScore_tt.VVV{1},'b', time, zScore_tt.VVN{1},'k',...
                            time, zScore_tt.VTN{1},'r', time, zScore_tt.VTV{1},'g');
                        hold on
                        sem_VVV = [zScore_tt.VVV{2}.', fliplr(zScore_tt.VVV{3}.')];
                        sem_VVN = [zScore_tt.VVN{2}.', fliplr(zScore_tt.VVN{3}.')];
                        sem_VTN = [zScore_tt.VTN{2}.', fliplr(zScore_tt.VTN{3}.')];
                        sem_VTV = [zScore_tt.VTV{2}.', fliplr(zScore_tt.VTV{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_VVV, 'b', time_bw, sem_VVN, 'k',...
                            time_bw, sem_VTN, 'r',time_bw, sem_VTV, 'g',...
                            'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Light Block', 'FontSize',16);
                        hold off

                        sgtitle(recSite, 'FontSize',20);
                        
                        % Save fig
                        zScoreFigPath = fullfile(save_path,[recSite,'_allTrials_BlockType.pdf']);
                        print(zScoreFigPath,'-dpdf','-painters','-loose');
                        
                    case 'stimulus'
                        subplot(1,2,1)
                        plot(time, zScore_tt.TTT{1}, 'b', time, zScore_tt.TTN{1},'k',...
                            time, zScore_tt.VTN{1},'r', time, zScore_tt.VTV{1},'g');
                        hold on
                        sem_TTT = [zScore_tt.TTT{2}.', fliplr(zScore_tt.TTT{3}.')];
                        sem_TTN = [zScore_tt.TTN{2}.', fliplr(zScore_tt.TTN{3}.')];
                        sem_VTN = [zScore_tt.VTN{2}.', fliplr(zScore_tt.VTN{3}.')];
                        sem_VTV = [zScore_tt.VTV{2}.', fliplr(zScore_tt.VTV{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_TTT, 'b', time_bw, sem_TTN, 'k',...
                            time_bw, sem_VTN, 'r',time_bw, sem_VTV, 'g',...
                            'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Tactile Stim', 'FontSize',16);
                        hold off

                        subplot(1,2,2)
                        plot(time, zScore_tt.VVV{1},'b', time, zScore_tt.VVN{1},'k',...
                            time,zScore_tt.TVN{1},'r', time, zScore_tt.TVT{1},'g');
                        hold on
                        sem_VVV = [zScore_tt.VVV{2}.', fliplr(zScore_tt.VVV{3}.')];
                        sem_VVN = [zScore_tt.VVN{2}.', fliplr(zScore_tt.VVN{3}.')];
                        sem_TVN = [zScore_tt.TVN{2}.', fliplr(zScore_tt.TVN{3}.')];
                        sem_TVT = [zScore_tt.TVT{2}.', fliplr(zScore_tt.TVT{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_VVV, 'b', time_bw, sem_VVN, 'k',...
                            time_bw, sem_TVN, 'r',time_bw, sem_TVT, 'g',...
                            'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Visual Stim', 'FontSize',16);
                        hold off

                        sgtitle(recSite, 'FontSize',20);
                        
                        % Save fig
                        zScoreFigPath = fullfile(save_path,[recSite,'_allTrials_StimType.pdf']);
                        print(zScoreFigPath,'-dpdf','-painters','-loose');
                end
            end
        
        end
        
        function plot_lickTrials (zScore_recSite, type, timeWindow, save_path)
            time = -1:0.01:2.49;
            for i = 1: size(zScore_recSite,1)
                zScore_tt = zScore_recSite{i,2};
                recSite = zScore_recSite{i,1};
                figure('Position', [0,0, 800, 500]); % set figure size
                switch type
                    case 'stimulus'
                        subplot(1,2,1)
        %                 plot(time, zScore_tt.TTT{1}, 'b', time, zScore_tt.TTV{1},'y',...
        %                 time, zScore_tt.VTV{1},'c', time, zScore_tt.VTT{1},'m');
                        plot(time, zScore_tt.TTT{1}, 'b',...
                        time, zScore_tt.VTV{1},'c', time, zScore_tt.VTT{1},'m');
                        hold on
                        sem_TTT = [zScore_tt.TTT{2}.', fliplr(zScore_tt.TTT{3}.')];
                        %sem_TTV = [zScore_tt.TTV{2}.', fliplr(zScore_tt.TTV{3}.')]; % exploration
                        sem_VTV = [zScore_tt.VTV{2}.', fliplr(zScore_tt.VTV{3}.')]; % compulsive licking 
                        sem_VTT = [zScore_tt.VTT{2}.', fliplr(zScore_tt.VTT{3}.')]; % rule error 
                        time_bw = [time, fliplr(time)];
        %                 fill(time_bw, sem_TTT, 'b', time_bw, sem_TTV, 'y',...
        %                 time_bw, sem_VTV, 'c',time_bw, sem_VTT, 'm',...
        %                 'FaceAlpha', 0.3, 'Linestyle', 'none');
                        fill(time_bw, sem_TTT, 'b',...
                        time_bw, sem_VTV, 'c',time_bw, sem_VTT, 'm',...
                        'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Tactile Stim', 'FontSize',16);
                        hold off

                        subplot(1,2,2)
        %                 plot(time, zScore_tt.VVV{1},'b', time, zScore_tt.VVT{1},'y',...
        %                 time,zScore_tt.TVT{1},'c', time, zScore_tt.TVV{1},'m');
                        plot(time, zScore_tt.VVV{1},'b',...
                        time,zScore_tt.TVT{1},'c', time, zScore_tt.TVV{1},'m');
                        hold on
                        sem_VVV = [zScore_tt.VVV{2}.', fliplr(zScore_tt.VVV{3}.')];
                        %sem_VVT = [zScore_tt.VVT{2}.', fliplr(zScore_tt.VVT{3}.')];
                        sem_TVT = [zScore_tt.TVT{2}.', fliplr(zScore_tt.TVT{3}.')];
                        sem_TVV = [zScore_tt.TVV{2}.', fliplr(zScore_tt.TVV{3}.')];
                        time_bw = [time, fliplr(time)];
        %                 fill(time_bw, sem_VVV, 'b', time_bw, sem_VVT, 'y',...
        %                 time_bw, sem_TVT, 'c',time_bw, sem_TVV, 'm',...
        %                 'FaceAlpha', 0.3, 'Linestyle', 'none');
                        fill(time_bw, sem_VVV, 'b',...
                        time_bw, sem_TVT, 'c',time_bw, sem_TVV, 'm',...
                        'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Visual Stim', 'FontSize',16);
                        hold off

                        sgtitle(recSite, 'FontSize',20);

                        % Save fig
                        zScoreFigPath = fullfile(save_path,[recSite,'_lickTrials_StimType_stimOnset.pdf']);
                        print(zScoreFigPath,'-dpdf','-painters','-loose');
                        
                    case 'lick direction'
                        subplot(1,2,1)
                        plot(time, zScore_tt.TTT{1}, 'b',...
                        time, zScore_tt.TVT{1},'c', time, zScore_tt.VTT{1},'m');
                        hold on
                        sem_TTT = [zScore_tt.TTT{2}.', fliplr(zScore_tt.TTT{3}.')];
                        sem_TVT = [zScore_tt.TVT{2}.', fliplr(zScore_tt.TVT{3}.')]; % compulsive licking 
                        sem_VTT = [zScore_tt.VTT{2}.', fliplr(zScore_tt.VTT{3}.')]; % rule error 
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_TTT, 'b',...
                        time_bw, sem_TVT, 'c',time_bw, sem_VTT, 'm',...
                        'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Right Lick', 'FontSize',16);
                        hold off

                        subplot(1,2,2)
                        plot(time, zScore_tt.VVV{1},'b',...
                        time,zScore_tt.VTV{1},'c', time, zScore_tt.TVV{1},'m');
                        hold on
                        sem_VVV = [zScore_tt.VVV{2}.', fliplr(zScore_tt.VVV{3}.')];
                        sem_VTV = [zScore_tt.VTV{2}.', fliplr(zScore_tt.VTV{3}.')];
                        sem_TVV = [zScore_tt.TVV{2}.', fliplr(zScore_tt.TVV{3}.')];
                        time_bw = [time, fliplr(time)];
                        fill(time_bw, sem_VVV, 'b',...
                        time_bw, sem_VTV, 'c',time_bw, sem_TVV, 'm',...
                        'FaceAlpha', 0.3, 'Linestyle', 'none');
                        xlim(timeWindow);
                        ylim([-0.5 3]);
                        xlabel('Time from stimulus onsets');
                        ylabel('Z-score');
                        title('Left Lick', 'FontSize',16);
                        hold off

                        sgtitle(recSite, 'FontSize',20);

                        % Save fig
                        zScoreFigPath = fullfile(save_path,[recSite,'_lickTrials_LickType_stimOnset.pdf']);
                        print(zScoreFigPath,'-dpdf','-painters','-loose');
                end        
            end
        end
        
        function plot_Heatmap(zScore_all, recSite, type, timeWindow, save_path)
            % Inputs
            %   zScore_all: zScore data for all recording sites
            %   recSite: a list of recording sites for plotting
            %   type: plot all,correct, or lick trials  
            %   time window (sec) ex: [-0.25 1] 0: stimulus onset
            %   save_path ex:
            %   'E:\CM_NeuralActivity_Analysis\zScore\GrandAvg'
            % Outputs
            %   heatmap figures for each recording site saved in a specific path
            
            % Load trial type index
            VariableNames = zScore_all.Properties.VariableNames;

            % select trial types
            switch type
                case 'all'
                    tt_numb = [4, 7, 6, 9, 11, 14, 12, 15]; %TTT TVT TTN TVN VTV VVV VTN VVN
                case 'correct'   
                    tt_numb = [4, 9, 12, 14]; %TTT TVN VTN VVV
                case 'lick'
                    tt_numb = [4 11 10 14 7 8]; %TTT VTV VTT VVV TVT TVV
            end
            % Convert time window to time bin window
            timeBinWindow = ((timeWindow+1)/0.01)+1;

            %zScore time window = -1 ~2.5
            for i = 1:length(recSite)
                HeatmapData = [];
                u = zScore_all(strcmp(zScore_all.recSite, recSite{i}), :);
                % Ctreate a new cell array for zScore heatmap
                for tt = 1:length(tt_numb)
                    v = [];
                    for ss = 1:height(u)
%                         isResponsive = u{ss, 16}{1,1};
%                         new_v = u{ss, tt_numb(tt)}{1,1}(:,isResponsive);
                        number_of_units = size(u{ss, 4}{1,1},2);
                        new_v = u{ss, tt_numb(tt)}{1,1};
                        if isempty(new_v)
                            % Apply NaN to miising data (because some sessions may not
                            % have all trial types)
                            new_v = NaN(350,number_of_units);
                        end                 
                        v = [v new_v];
                    end
                    HeatmapData = [HeatmapData; v];
                end

                HeatmapData = HeatmapData.';
                tHitScore = [];
                for j = 1: size(HeatmapData,1)
                    new_score = mean(HeatmapData(j,101:150)); %sorted by tHit 0~0.5 s
                    tHitScore = [tHitScore; new_score];
                end
                HeatmapData = [HeatmapData tHitScore];
                sorted_data = sortrows(HeatmapData, size(HeatmapData,2), 'descend');

                % Plot
                xTickLables = [repmat(' ',99,1); '0'; repmat(' ',99,1);'1'; repmat(' ',99,1); '2'; repmat(' ',50,1)]; 
                yTickLables = [{num2str(size(sorted_data,1))}; repmat({' '},size(sorted_data,1)-2,1); {'1'}]; 

                switch type
                    case 'all'
                        figure('Position', [0,0, 800,1000]);clf % set figure size
                        for tt = 1:length(tt_numb)  
                            subplot(2,4,tt)
                            if tt == 1 || tt == 5
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)', 'YLabel','Unit', 'ColorbarVisible', 'off','GridVisible','off',...
                                'Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',yTickLables);
                            elseif tt == 4 || tt == 8
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6],'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)','GridVisible','off','Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',repmat(' ',size(sorted_data,1),1));
                            else 
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6],'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)', 'ColorbarVisible', 'off','GridVisible','off','Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',repmat(' ',size(sorted_data,1),1));
                            end
                        end
                        sgtitle(recSite{i},'FontSize',10);
                        fig = gcf;
                        HMFigPath = fullfile(save_path,[recSite{i},'_allTrials_HM.pdf']);
                        print(HMFigPath,'-dpdf','-painters','-loose');

                    case 'correct'
                        figure('Position', [0,0, 800,500]);clf % set figure size
                        for tt = 1:length(tt_numb)  
                            subplot(1,4,tt)
                            if tt == 1
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)', 'YLabel','Unit', 'ColorbarVisible', 'off','GridVisible','off',...
                                'Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',yTickLables);
                            elseif tt == 4 
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)','GridVisible','off','Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',repmat(' ',size(sorted_data,1),1));
                            else 
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)', 'ColorbarVisible', 'off','GridVisible','off','Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',repmat(' ',size(sorted_data,1),1));
                            end
                        end
                        sgtitle(recSite{i},'FontSize',10);
                        fig = gcf;
                        HMFigPath = fullfile(save_path,[recSite{i},'_correctTrials_HM.pdf']);
                        print(HMFigPath,'-dpdf','-painters','-loose');

                    case 'lick'
                        figure('Position', [0,0, 800,1000]);clf % set figure size
                        for tt = 1:length(tt_numb)  
                            subplot(2,3,tt)
                            if tt == 1 || tt == 4
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)', 'YLabel','Unit', 'ColorbarVisible', 'off','GridVisible','off',...
                                'Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',yTickLables);
                            elseif tt == 3 || tt == 6
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)','GridVisible','off','Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',repmat(' ',size(sorted_data,1),1));
                            else 
                                h = heatmap(sorted_data(:, 1+(tt-1)*350 :tt*350),'Colormap',jet, 'ColorLimits',[-4 6], 'Xlimits', timeBinWindow,...
                                'XLabel','Time (s)', 'ColorbarVisible', 'off','GridVisible','off','Title', VariableNames{tt_numb(tt)}, ...
                                'XDisplayLabels',xTickLables,'YDisplayLabels',repmat(' ',size(sorted_data,1),1));
                            end
                        end
                        sgtitle(recSite{i},'FontSize',10);
                        fig = gcf;
                        HMFigPath = fullfile(save_path,[recSite{i},'_lickTrials_HM.pdf']);
                        print(HMFigPath,'-dpdf','-painters','-loose');                       
                end 
            end     
        end
        
    end
end
