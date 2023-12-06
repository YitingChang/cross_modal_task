%% Create a class
classdef tCorrectTrials
    % functions for analyzing tactile correct trials: tHit and tCR
    
    methods(Static)
        function GetAmplitude(MouseName, seDir)
            % Inputs
            %   MouseName    Mouse name
            %   seDir        Directory of session explore
            % Outputs
            %   Normalized responses of tHit and tCR trials for each units
            %   Analysis window: 0~ 0.25 second
            %   Baseline window: -1~0 second
            %   Normalization: (response - mean(baseline))/mean(baseline)
            
            seFilePaths = MBrowse.Files(seDir);

            tCorrectTrials_summary = [];
            for i = 1: length(seFilePaths)
                load(seFilePaths{i})
                tCorrectTrials = [];
                spkRs_baseline = CM.Analysis.GetSpikeRate(se, 1, -1, 0); % baseline = -1~0 sec
                spkRs_response = CM.Analysis.GetSpikeRate(se, 0.25, 0, 0.25); % response = 0 ~ 0.25 sec
                trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                     'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
                ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
                'VariableNames', trialTypes);
                isResponsive = CM.Analysis.GetResponsiveUnit(se, [-1, 2.5]); 
                
                 % select unit
                for u = 2:size(spkRs_response, 2) % first column is time and unit starts from 2nd column
                    if isResponsive(u-1) 
                        sb = cellfun(@(x) x/4, spkRs_baseline{:,u}); % mean baseline
                        sr = cell2mat(spkRs_response{:,u});
                        % select trial types
                        sr_TTT = sr(ttInd.TTT);
                        sr_VTN = sr(ttInd.VTN);
                        normalized_TTT = (mean(sr_TTT)- mean(sb))/mean(sb);
                        normalized_VTN = (mean(sr_VTN)- mean(sb))/mean(sb);
                        isDiff_TTT = ttest2(sb,sr_TTT);
                        isDiff_VTN = ttest2(sb,sr_VTN);

                        % summary
                        tCorrectTrials{1, u-1} = normalized_TTT;
                        tCorrectTrials{2, u-1} = normalized_VTN;
                        tCorrectTrials{3, u-1} = isDiff_TTT;
                        tCorrectTrials{4, u-1} = isDiff_VTN;
                    else
                        % Stop if it is not responsive
                        continue
                    end

                end    
                
                if ~isempty(tCorrectTrials)
                    % Load session information 
                    MouseName = se.userData.sessionInfo.MouseName;
                    seshDate = se.userData.sessionInfo.seshDate;
                    recSite = se.userData.sessionInfo.recSite;

                    % 
                    new_tCorrectTrials = [{MouseName} {seshDate} {recSite} {tCorrectTrials}];

                    % Concatenation 
                    tCorrectTrials_summary = [tCorrectTrials_summary; new_tCorrectTrials];
                else
                    % Stop if there is no responsive unit in this session
                    continue
                end
                clearvars -except seFilePaths tCorrectTrials_summary MouseName
            end

            % Convert to a table
            VarNames = {'MouseName', 'seshDate', 'recSite', 'tCorrectTrials_summary'};
            tCorrectTrials_summary = cell2table(tCorrectTrials_summary,'VariableNames', VarNames);
            
            % Save tCorrectTrials_summary summary for each mouse 
            tCorrectPaths = fullfile('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Amplitude',[MouseName,'_tCorrectTrials']);
            save(tCorrectPaths, 'tCorrectTrials_summary');
        end
        
        function allConcat(tCorrectAmpDir)
            tCorrectAmpPaths = MBrowse.Files(tCorrectAmpDir);
            tCorrectAmp_all = table();
            for i = 1: length(tCorrectAmpPaths)
                load(tCorrectAmpPaths{i})
                tCorrectAmp_all = [tCorrectAmp_all; tCorrectTrials_summary];
            end
                tCorrectAmpSavePath = fullfile(tCorrectAmpDir, 'tCorrectAmp');
                save(tCorrectAmpSavePath, 'tCorrectAmp_all');
        end
        
        function ranSelect_amp(tCorrectAmp_tetrode, n)
            % Randomly select n number of units from tetrode dataset 

            % Concatenation
            tCorrectAmp = [];
            for j = 1:size(tCorrectAmp_tetrode,1) % loop each session
                tCorrectAmp = [tCorrectAmp tCorrectAmp_tetrode{j,4}{1,1}];
            end
            rng(19); % control random number generation
            randomInd = randperm(size(tCorrectAmp,2),n); % randomly select n number from the tototal number of units

            tCorrectAmp_selected = tCorrectAmp(:,randomInd);

            recSite = tCorrectAmp_tetrode.recSite{1};
            tCorrectAmp_tetrode_selected = [{'Random'} {'Random'} {recSite} {tCorrectAmp_selected} ]; 
            % Convert to a table
            VarNames = {'MouseName', 'seshDate', 'recSite', 'tCorrectTrials_summary'};
            tCorrectAmp_ranSelect = cell2table(tCorrectAmp_tetrode_selected,'VariableNames', VarNames);

            % Save btInd summary for each mouse 
            tCorrectAmpPaths = fullfile('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Amplitude',...
                ['Random', '_', recSite, '_tCorrectAmp']);
            save(tCorrectAmpPaths, 'tCorrectAmp_ranSelect');            
        end
        
        function tetrode2siConcat_amp(tCorrectAmpDir, tCorrectAmp_SiProbe)
            % Input
            %   btIndDir  where to find block type index of tetrode dataset
            %   btInd_SiProbe block type index of silicone probe dataset
            % Concatenation 
            tCorrectAmpPaths = MBrowse.Files(tCorrectAmpDir);
            tCorrectAmp_all = tCorrectAmp_SiProbe;
            for i = 1: length(tCorrectAmpPaths)
                load(tCorrectAmpPaths{i})
                tCorrectAmp_all = [tCorrectAmp_all; tCorrectAmp_ranSelect];
            end

            % Save tCorrectAmp summary for all mice 
            save('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Amplitude\tCorrectAmp_mixed', 'tCorrectAmp_all');
        end
        
        function tCorrectAmp_recSite = GetAmpSummary_recSite(tCorrectTrials_all, recSite)
%             Inputs
%                 tHit and tCR of units from region of interests (recording sites)
%             Outputs 
%                 1) tHit and tCR for different significant groups (sig_tHit, sig_tCR, sig_tHitCR, nonsig_tHitCR)
%                 2) The percentage of four differnet significant groups
%                 3) x and y label names for plotting
            
            for i = 1: length(recSite)
                t1 = tCorrectTrials_all(strcmp(tCorrectTrials_all.recSite, recSite{i}), :); % select recording sites
                tHit_sig_tHit = [];
                tCR_sig_tHit = [];
                tHit_sig_tCR = [];
                tCR_sig_tCR = [];
                tHit_sig_tHitCR = [];
                tCR_sig_tHitCR = [];
                tHit_nonsig_tHitCR = [];
                tCR_nonsig_tHitCR= [];
                for j = 1:size(t1,1) % loop each session
                    t2 = cell2mat(t1{j,4}{1,1}); % get tCorrect summary

                    % remove unit having NaN and inf
                    t2(:,any(isnan(t2), 1)) = [];
                    t2(:,any(isinf(t2), 1)) = [];

                    % tHit- blue, tCR- red, both - purple
                    % Get all units
                    tHit_all = t2(1,:);
                    tCR_all = t2(2,:);
                    % Get units with significant tHit
                    new_tHit_sig_tHit = tHit_all(logical(t2(3,:))&~logical(t2(4,:)));
                    new_tCR_sig_tHit = tCR_all(logical(t2(3,:))&~logical(t2(4,:)));
                    % Get units with significant tCR
                    new_tHit_sig_tCR = tHit_all(~logical(t2(3,:))&logical(t2(4,:)));
                    new_tCR_sig_tCR = tCR_all(~logical(t2(3,:))&logical(t2(4,:)));        
                    % Get units with significant tHit & tCR
                    new_tHit_sig_tHitCR = tHit_all(logical(t2(3,:))&logical(t2(4,:)));
                    new_tCR_sig_tHitCR = tCR_all(logical(t2(3,:))&logical(t2(4,:)));  
                    % Get units with non-significant tHit & tCR
                    new_tHit_nonsig_tHitCR = tHit_all(~logical(t2(3,:))&~logical(t2(4,:)));
                    new_tCR_nonsig_tHitCR = tCR_all(~logical(t2(3,:))&~logical(t2(4,:)));  

                    tHit_sig_tHit = [tHit_sig_tHit new_tHit_sig_tHit];
                    tCR_sig_tHit = [tCR_sig_tHit new_tCR_sig_tHit];
                    tHit_sig_tCR = [tHit_sig_tCR new_tHit_sig_tCR];
                    tCR_sig_tCR = [tCR_sig_tCR new_tCR_sig_tCR];
                    tHit_sig_tHitCR = [tHit_sig_tHitCR new_tHit_sig_tHitCR];
                    tCR_sig_tHitCR = [tCR_sig_tHitCR new_tCR_sig_tHitCR];
                    tHit_nonsig_tHitCR = [tHit_nonsig_tHitCR new_tHit_nonsig_tHitCR];
                    tCR_nonsig_tHitCR= [tCR_nonsig_tHitCR new_tCR_nonsig_tHitCR];
                end
                
                % percentage of significant units in recording sites
                total_units = length(tHit_sig_tHit)+ length(tHit_sig_tHitCR)+length(tHit_sig_tCR)+length(tHit_nonsig_tHitCR);
                perc_sig_tHit = (length(tHit_sig_tHit)+ length(tHit_sig_tHitCR))/ total_units;
                perc_sig_tCR = (length(tHit_sig_tCR)+ length(tHit_sig_tHitCR))/ total_units; 
                sig_summary = {[length(tHit_sig_tHitCR),length(tHit_sig_tCR);length(tHit_sig_tHit),length(tHit_nonsig_tHitCR)]./total_units}; 
                xlabelNames = {strcat('tHit (', num2str(round(perc_sig_tHit,2)*100), '%)')};
                ylabelNames = {strcat('tCR (', num2str(round(perc_sig_tCR,2)*100), '%)')};
                
                tCorrect_summary = [{tHit_sig_tHit} {tCR_sig_tHit} {tHit_sig_tCR} {tCR_sig_tCR}...
                                    {tHit_sig_tHitCR} {tCR_sig_tHitCR} {tHit_nonsig_tHitCR} {tCR_nonsig_tHitCR}...
                                    {sig_summary} {xlabelNames} {ylabelNames}];
                                
                tCorrectAmp_recSite{i,1} = recSite{i};
                tCorrectAmp_recSite{i,2} = tCorrect_summary; 
            end
        end
        
        function plot_Amp (tCorrectAmp_recSite, recSite, save_path)
                figure('Position', [0,0, 600, 300*round(length(recSite)/2)]); % set figure size
                for i = 1: length(recSite)
                    subplot(round(length(recSite)/2),2,i)
                    hold on
                    % different groups are colored 
%                     scatter(tCorrectAmp_recSite{i,2}{1},tCorrectAmp_recSite{i,2}{2},10, 'b', 'filled') % significant tHit
%                     scatter(tCorrectAmp_recSite{i,2}{3},tCorrectAmp_recSite{i,2}{4},10, 'r', 'filled') % significant tCR
%                     scatter(tCorrectAmp_recSite{i,2}{5},tCorrectAmp_recSite{i,2}{6},10, 'g', 'filled') % significant tHit and tCR
%                     scatter(tCorrectAmp_recSite{i,2}{7},tCorrectAmp_recSite{i,2}{8},10, 'k') % non significant tHit and tCR
%                     xline(0,'--k');
%                     yline(0, '--k');
%                     line([-2.5 10],[-2.5 10],'Color','black','LineStyle','--')
%                     xlim([-2 10]);
%                     ylim([-2 10]);
%                     xticks([-2:2:10]);
%                     yticks([-2:2:10]);
%                     text(1.5,9, 'Hs          Hn','FontSize',6)
%                     text(0.5,7.5, ['Cs';'Cn'],'FontSize',6)
%                     text(1.5,7.5,num2str(round(tCorrectAmp_recSite{i,2}{9}{1,1},2)),'FontSize',6)
%                     xlabel(tCorrectAmp_recSite{i,2}{10})
%                     ylabel(tCorrectAmp_recSite{i,2}{11})
%                     title(recSite{i})
%                     hold off
%                 % Save fig
%                 tCorrectAmpFigPath = fullfile(save_path,'tCorrectAmp.pdf');
%                 print(tCorrectAmpFigPath,'-dpdf','-painters','-loose')
                    % non-colored
                    scatter(tCorrectAmp_recSite{i,2}{1},tCorrectAmp_recSite{i,2}{2},10, 'k') % significant tHit
                    scatter(tCorrectAmp_recSite{i,2}{3},tCorrectAmp_recSite{i,2}{4},10, 'k') % significant tCR
                    scatter(tCorrectAmp_recSite{i,2}{5},tCorrectAmp_recSite{i,2}{6},10, 'k') % significant tHit and tCR
                    scatter(tCorrectAmp_recSite{i,2}{7},tCorrectAmp_recSite{i,2}{8},10, 'k') % non significant tHit and tCR
                    xline(0,'--b');
                    yline(0, '--b');
                    line([-2.5 10],[-2.5 10],'Color','blue','LineStyle','--')
                    xlim([-2 10]);
                    ylim([-2 10]);
                    xticks([-2:2:10]);
                    yticks([-2:2:10]);
                    xlabel('tHit z-score')
                    ylabel('tCR z-score')
                    title(recSite{i})
                    hold off
                    % Save fig
                    tCorrectAmpFigPath = fullfile(save_path,'tCorrectAmp_non-colored.pdf');
                    print(tCorrectAmpFigPath,'-dpdf','-painters','-loose')
                end

        end
        
        function AmpStatistics = GetAmpStatistics_recSite(tCorrectAmp_all, recSite)
            
            for i = 1: length(recSite)
                t1 = tCorrectAmp_all(strcmp(tCorrectAmp_all.recSite, recSite{i}), :); % select recording sites
                tHit = [];
                tCR = [];
                for j = 1:size(t1,1) % loop each session
                    t2 = cell2mat(t1{j,4}{1,1}); % get tCorrect summary

                    % remove unit having NaN and inf
                    t2(:,any(isnan(t2), 1)) = [];
                    t2(:,any(isinf(t2), 1)) = [];

                    % Get all units
                    tHit_all = t2(1,:);
                    tCR_all = t2(2,:);

                    tHit = [tHit tHit_all];
                    tCR = [tCR tCR_all];
                end

                % Using the bootstrap to test whether two sets of data (tHit vs. tCR)
                % are dufferent in their means
                unitNum = length(tHit);
                % Aggrgate the data, draw bootstrap samples, and compute the resulting
                % differences in means
                rng(53); 
                dist = bootstrp(10000, @(x) mean(x(1:unitNum))- mean(x(unitNum+1:2*unitNum)), [tHit tCR]); 

                % Compute the actual observed difference in means as well as the
                % corresponding p-value.
                val = mean(tHit) - mean(tCR);
                pval = sum(abs(dist)>abs(val))/length(dist);

                AmpStatistics{i,1} = recSite{i};
                AmpStatistics{i,2} = unitNum;
                AmpStatistics{i,3} = val;
                AmpStatistics{i,4} = pval;
            end
            % Convert to a table
            VarNames = {'recSite' 'n' 'diffMean' 'pValue'};
            AmpStatistics = cell2table(AmpStatistics,'VariableNames', VarNames);
        end
        
        function GetUnitSummary(MouseName, seDir)
            % Input
            % spike data (spkData) from MSessionExplore
       
            % Output
            % MouseName, seshDate, recSite, zScore
            % zScore: 
            % rows - mean, lower bound and higher bound of 95%CI of each time stamp for
            % each unit (bin: 0.025 second, timeWindow: 0~2 second)
            % columns - tHit (TTT), tCR(VTN)
            
            seFilePaths = MBrowse.Files(seDir);
            zScore = [];
            for i = 1: length(seFilePaths)
                load(seFilePaths{i})
                
                timeWindow = [0 2];
                bin = 0.025;
                spkRs = CM.Analysis.GetSpikeRate(se, bin, timeWindow(1), timeWindow(2)); 
                normalized_spkRs = CM.Analysis.zScore(se, bin, spkRs);
                isResponsive = CM.Analysis.GetResponsiveUnit(se, [-1 2.5]); 
                normalized_responsive_spkRs = normalized_spkRs(:,isResponsive); % select responsive units
                trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
                ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
                'VariableNames', trialTypes);
                % select trial types: TTT (tHit) and VTN (tCR)
                tt_ind = [1 9]; % TTT: 1st column, VTN: 9th column
    
                sig_onset_session = NaN(2, size(normalized_responsive_spkRs, 2));
                
                for tt = 1:size(tt_ind, 2)
                    ind = ttInd{:, tt_ind(tt)};
                    mean_tt = [];
                    lbCI_tt = [];
                    hbCI_tt = [];
                    sig_ind_tt = [];
                    parfor unit = 1 : size(normalized_responsive_spkRs, 2) % parallel for loop 
                        s = normalized_responsive_spkRs{:,unit}; % Select unit
                        s = s(ind); % Select trial types
                        s2 = reshape(cell2mat(s), diff(timeWindow)/bin, size(s,1));
                        s2 = num2cell(s2,2);
                        %%%% bootstrapping %%%%
                        rng(22); % control random number generation
                        nboot = 1000; % nboot
                        bootstat = cellfun(@(x) bootstrp(nboot,@mean,x), s2, 'UniformOutput', false);
                        bootstat_sorted = cellfun(@(x) sort(x), bootstat, 'UniformOutput', false);
                        mean_tt_unit =  cellfun(@(x) mean(x), bootstat_sorted,'UniformOutput', false);
                        lbCI_tt_unit = cellfun(@(x)x(nboot*0.025), bootstat_sorted,'UniformOutput', false);
                        hbCI_tt_unit = cellfun(@(x)x(nboot*0.975), bootstat_sorted,'UniformOutput', false);
                        sig_ind_tt_unit = cell2mat(lbCI_tt_unit).*cell2mat(hbCI_tt_unit) >0; % if 95%CI contains 0, lbCI*hbCi will be negative                          

                        % Add each unit trace to a table 
                        mean_tt = [mean_tt, mean_tt_unit];
                        lbCI_tt = [lbCI_tt, lbCI_tt_unit];
                        hbCI_tt = [hbCI_tt, hbCI_tt_unit];
                        sig_ind_tt = [sig_ind_tt, sig_ind_tt_unit]; 
                    end 
                    spkRate_tt{1, tt} = mean_tt;
                    spkRate_tt{2, tt} = lbCI_tt;
                    spkRate_tt{3, tt} = hbCI_tt;
                    spkRate_tt{4, tt} = sig_ind_tt;
                end
                
                % Load session information 
                MouseName = se.userData.sessionInfo.MouseName;
                seshDate = se.userData.sessionInfo.seshDate;
                recSite = se.userData.sessionInfo.recSite;

                % Combine with zScore
                new_zScore = [{MouseName} {seshDate} {recSite} {spkRate_tt}];

                % Concatenation 
                zScore = [zScore; new_zScore];
                clearvars -except seFilePaths zScore MouseName bin timeWindow          
            end
            % Convert to a table
            VarNames = {'MouseName' 'seshDate' 'recSite' 'zScore'};
            zScore = cell2table(zScore,'VariableNames', VarNames);
            % Save zScore
            zScorePaths = fullfile('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Onset\zScore', ...
                        [MouseName,'_zScore']);
            save(zScorePaths, 'zScore');            
        end  
        
        function GetOnset(zScoreDir)
            % Inputs
            %   MouseName    Mouse name
            %   zScoreDir    Directory of zScore for getting the onsets of response 
            % Outputs
            %   Get response onsets for each unit
            %   Index of significance: For each time bin, whether 95%CI of activity
            %   covers 0. If not, it is significant for that time bin.
            %   The Onset of response is determined by the first time bin of three consecutive
            %   significant time bins after stimulus onset.

            zScoreFilePaths = MBrowse.Files(zScoreDir);
            time = 0:0.025:2-0.025;
            tCorrectOnset_all = [];
            for i = 1: length(zScoreFilePaths)
                load(zScoreFilePaths{i})
                sig_onset = [];
                for session = 1: size(zScore,1)
                    zScore_session = zScore{session, 4}{1,1}; 
                    sig_onset_session = NaN(2, size(zScore_session{4,1}, 2));
                    for tt = 1:size(zScore_session, 2)
                        sig = zScore_session{4,tt}; % significant index is 4th row
                        for unit = 1:size(sig, 2)
                            sig_unit = sig(:,unit);
                            sig_sum = sum([sig_unit, [sig_unit(2:end);0],[sig_unit(3:end);0;0]], 2);
                            sig_bin = find(sig_sum==3, 1); 
                            if ~isempty(sig_bin)
                                sig_onset_session(tt,unit) = time(sig_bin);
                            end
                        end
                    end
                    % Load session information 
                    MouseName = zScore{session, 1}{1,1};
                    seshDate = zScore{session, 2}{1,1};
                    recSite = zScore{session, 3}{1,1};
                    % Session summary
                    new_sig_onset = [{MouseName} {seshDate} {recSite} {sig_onset_session}];
                    % Concatenation 
                    sig_onset = [sig_onset; new_sig_onset];
                end
                tCorrectOnset_all = [tCorrectOnset_all; sig_onset];
                clearvars -except zScoreFilePaths tCorrectOnset_all time
            end
            % Convert to a table
            VarNames = {'MouseName' 'seshDate' 'recSite' 'sig_onset'};
            tCorrectOnset_all = cell2table(tCorrectOnset_all,'VariableNames', VarNames);

            % Save onset summary for all mouse 
            onsetPaths = fullfile('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Onset', 'tCorrectOnset_all');
            save(onsetPaths, 'tCorrectOnset_all');     
        end   
        
        function ranSelect_onset(tCorrectOnset_tetrode, n)
            % Randomly select n number of units from tetrode dataset 
            % Concatenation
            tCorrectOnset = [];
            for j = 1:size(tCorrectOnset_tetrode,1) % loop each session
                tCorrectOnset = [tCorrectOnset tCorrectOnset_tetrode{j,4}{1,1}];
            end
            rng(19); % control random number generation
            randomInd = randperm(size(tCorrectOnset,2),n); % randomly select n number from the tototal number of units

            tCorrectOnset_selected = tCorrectOnset(:,randomInd);

            recSite = tCorrectOnset_tetrode.recSite{1};
            tCorrectOnset_tetrode_selected = [{'Random'} {'Random'} {recSite} {tCorrectOnset_selected} ]; 
            % Convert to a table
            VarNames = {'MouseName', 'seshDate', 'recSite', 'sig_onset'};
            tCorrectOnset_ranSelect = cell2table(tCorrectOnset_tetrode_selected,'VariableNames', VarNames);

            % Save btInd summary for each mouse 
            tCorrectOnsetPaths = fullfile('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Onset',...
                ['Random', '_', recSite, '_tCorrectOnset']);
            save(tCorrectOnsetPaths, 'tCorrectOnset_ranSelect');            
        end
        
        function tetrode2siConcat_onset(tCorrectOnsetDir, tCorrectOnset_SiProbe)
            % Input
            %   tCorrectOnsetDir:  where to find onsets of tetrode dataset
            %   tCorrectOnset_SiProbe: block type index of silicone probe dataset
            % Concatenation 
            tCorrectOnsetPaths = MBrowse.Files(tCorrectOnsetDir);
            tCorrectOnset_all = tCorrectOnset_SiProbe;
            for i = 1: length(tCorrectOnsetPaths)
                load(tCorrectOnsetPaths{i})
                tCorrectOnset_all = [tCorrectOnset_all; tCorrectOnset_ranSelect];
            end

            % Save btInx summary for all mice 
            save('E:\CM_NeuralActivity_Analysis\tCorrectTrials\Onset\tCorrectOnset_mixed', 'tCorrectOnset_all');
        end
        
        
        function tCorrectOnset_recSite = GetOnsetSummary_recSite(onset_all, recSite)
            for i = 1: length(recSite)
                t1 = onset_all(strcmp(onset_all.recSite, recSite{i}), :); % select recording sites
                t3 = [];
                for j = 1:size(t1,1) % loop each session
                    t2 = t1{j,4}{1,1}; % get onset summary
                    % remove unit having NaN and inf
                    t2(:,any(isnan(t2), 1)) = [];
                    t2(:,any(isinf(t2), 1)) = [];
                    t3 = [t3 t2];
                end
                tCorrectOnset_recSite{i, 1} = recSite{i};
                tCorrectOnset_recSite{i, 2} = [t3];
            end
        end
        
        function OnsetStatistics = GetOnsetStatistics_recSite(onset_all, recSite)
            for i = 1: length(recSite)
                t1 = onset_all(strcmp(onset_all.recSite, recSite{i}), :); % select recording sites
                t3 = [];
                for j = 1:size(t1,1) % loop each session
                    t2 = t1{j,4}{1,1}; % get onset summary
                    % remove unit having NaN and inf
                    t2(:,any(isnan(t2), 1)) = [];
                    t2(:,any(isinf(t2), 1)) = [];
                    t3 = [t3 t2];
                end

                % Using the bootstrap to test whether two sets of data (tHit vs. tCR onsets)
                % are dufferent in their means
                tHit = t3(1,:);
                tCR = t3(2,:);
                unitNum = length(tHit);
                % Aggrgate the data, draw bootstrap samples, and compute the resulting
                % differences in means
                rng(48); 
                dist = bootstrp(10000, @(x) mean(x(1:unitNum))- mean(x(unitNum+1:2*unitNum)), [tHit tCR]); 

                % Compute the actual observed difference in means as well as the
                % corresponding p-value.
                val = mean(tHit) - mean(tCR);
                pval = sum(abs(dist)>abs(val))/length(dist);

                OnsetStatistics{i,1} = recSite{i};
                OnsetStatistics{i,2} = unitNum;
                OnsetStatistics{i,3} = val;
                OnsetStatistics{i,4} = pval; 
            end
            % Convert to a table
            VarNames = {'recSite' 'n' 'diffMean' 'pValue'};
            OnsetStatistics = cell2table(OnsetStatistics,'VariableNames', VarNames);
        end

        
        function plot_Onset(tCorrectOnset_recSite, recSite, save_path)
            figure('Position', [0,0, 400, 200*round(length(recSite)/2)]); % set figure size
            for i = 1: length(recSite)
                tHit = tCorrectOnset_recSite{i,2}(1,:);
                tCR = tCorrectOnset_recSite{i,2}(2,:);
                subplot(round(length(recSite)/2),2,i)
                hold on
                scatter(tHit,tCR,10, 'k', 'filled')
                line([-2.5 10],[-2.5 10],'Color','black','LineStyle','--')
                xlim([0 2]);
                ylim([0 2]);
                xlabel('Time from stimulus onset (s) tHit')
                ylabel('Time from stimulus onset (s) tCR')
                xticks([0:0.4:2]);
                yticks([0:0.4:2]);
                title(recSite{i})
                hold off
            end
            % Save fig
            tCorrectOnsetFigPath = fullfile(save_path,'tCorrectOnset_tetrode.pdf');
            print(tCorrectOnsetFigPath,'-dpdf','-painters','-loose')
        end
        
    end
            
        
end


