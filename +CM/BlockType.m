%% Analysis of block type index
%  Create a class
classdef BlockType
    % functions for analyzing block type index
    
    methods(Static)
        function GetBTind(MouseName, seDir)
            % Inputs
            %   MouseName    Mouse name
            %   seDir        Directory of session explore
            % Output
            %   Sessions at least 4 blocks included
            %   The sumary of block type index includes mean, lower ci, upper ci, and significant
            %   results for each unit - Bootstrapping
            %   The drifting_summary: 1 is drifting and 0 is not drifting- ANOVA
            
            seFilePaths = MBrowse.Files(seDir);
            btInds = [];
            for i = 1: length(seFilePaths)
                load(seFilePaths{i})
                % Quality control- session
                [tBlock_onset,vBlock_onset,tBlock_endLength,vBlock_endLength] = CM.QualityControl.GetBlockInfo(se);  
%                 %%%%%%%% Create index for block numbers %%%%%%%%
%                 trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
%                              'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
%                 ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
%                         'VariableNames', trialTypes);
% 
%                 BlockChangeInd = ismember(se.GetColumn('behavValue', 'blockType'), 'Change') ;
% 
%                 tBlockInd = sum(ttInd{:,1:6},2); 
%                 vBlockInd = sum(ttInd{:,7:12},2);
%                 % Assign block type to the cahnge trials  
%                 for t = 2: length(tBlockInd)
%                     if BlockChangeInd(t) == 1 && tBlockInd(t-1) ==1
%                         tBlockInd(t) = 1;
%                     elseif BlockChangeInd(t) == 1 && vBlockInd(t-1) ==1
%                         vBlockInd(t) = 1;
%                     end 
%                 end
% 
%                 % Find block onset 
%                 [tBlock_onset,~] = find(diff([0;tBlockInd])==1); % onset = 1; offset = -1
%                 [vBlock_onset,~] = find(diff([0;vBlockInd])==1);
% 
%                 % Get how many trials in the last tblock and vBlock
%                 tBlock_endLength = length(tBlockInd) - tBlock_onset(length(tBlock_onset));
%                 vBlock_endLength = length(vBlockInd) - vBlock_onset(length(vBlock_onset));

                % Stop if number of blocks are less than 2 
                if length(tBlock_onset) <2 || length(vBlock_onset) <2 
                    continue
                elseif length(tBlock_onset) == 2 && length(vBlock_onset) == 2 
                    % Stop if it is 4 blocks but the last block is shorter than 30 trials 
                    if tBlock_endLength < 30 || vBlock_endLength < 30
                        continue
                    end
                end

                % Assign block numbers within each block type (1st tBlock =1, 2nd tBlock =2...)
                tBlockInd_number = tBlockInd;
                vBlockInd_number = vBlockInd;
                for t = 1: length(tBlockInd_number) % each trial
                    for n = 2: length(tBlock_onset) % each tBlock
                        if t >= tBlock_onset(n) && tBlockInd_number(t) == n-1
                            tBlockInd_number(t) = tBlockInd_number(t) + 1;
                        end
                    end   
                    for m = 2: length(vBlock_onset)
                        if t >= vBlock_onset(m) && vBlockInd_number(t) == m-1
                            vBlockInd_number(t) = vBlockInd_number(t) + 1;
                        end
                    end
                end

                %%%%%%%% Get drifting and block index for each unit %%%%%%%%
                spkRs = CM.Analysis.GetSpikeRate(se, 1, -1, 0); % (Hz)
                baseLine = table2cell(spkRs(:,2:end));
                tBlockInds = logical(repmat(tBlockInd,1,size(baseLine,2)));
                vBlockInds = logical(repmat(vBlockInd,1,size(baseLine,2)));

                tBaseLine = reshape(baseLine(tBlockInds),[sum(tBlockInd), size(baseLine,2)]);
                vBaseLine = reshape(baseLine(vBlockInds),[sum(vBlockInd), size(baseLine,2)]);


                for n = 1: length(tBlock_onset)
                    tBlockInds_block = repmat(ismember(tBlockInd_number,[n]), 1, size(baseLine,2));
                    tBaseLine_block(n) = {reshape(baseLine(tBlockInds_block),[sum(ismember(tBlockInd_number,[n])), size(baseLine,2)])};
                end 
                for m = 1: length(vBlock_onset)
                    vBlockInds_block = repmat(ismember(vBlockInd_number,[m]), 1, size(baseLine,2));
                    vBaseLine_block(m) = {reshape(baseLine(vBlockInds_block),[sum(ismember(vBlockInd_number,[m])), size(baseLine,2)])};
                end 


                for j = 1:size(baseLine, 2)  % each unit
                    tBlock_unit = [];
                    vBlock_unit = [];
                    tBlock_grp = [];
                    vBlock_grp = [];

                    % ANOVA: compare baseline activity among different tblocks 
                    for n = 1:length(tBlock_onset)
                        tBlock_unit = [tBlock_unit; tBaseLine_block{1,n}(:,j)]; % Concatenation
                        tBlock_grp = [tBlock_grp; repmat(n,size(tBaseLine_block{1,n},1) ,1)]; % Create group number
                    end
                    % ANOVA: compare baseline activity among different vblocks 
                    for m = 1:length(vBlock_onset)
                        vBlock_unit = [vBlock_unit; vBaseLine_block{1,m}(:,j)]; % Concatenation
                        vBlock_grp = [vBlock_grp; repmat(m,size(vBaseLine_block{1,m},1) ,1)]; % Create group number
                    end

                    p1 = anova1(cell2mat(tBlock_unit), tBlock_grp, 'off');
                    p2 = anova1(cell2mat(vBlock_unit), vBlock_grp, 'off');

                    if p1 > 0.05 && p2 > 0.05
                        drifting_summary{1,j} = 0;
                    else
                        drifting_summary{1,j} = 1;
                    end

                    % Index for block type
                    n = 1000; % nboot
                    t_bootstat = bootstrp(n,@mean,cell2mat(tBaseLine(:,j)));
                    v_bootstat = bootstrp(n,@mean,cell2mat(vBaseLine(:,j)));
                    btInd = (t_bootstat - v_bootstat)./(t_bootstat + v_bootstat); % Index for block type
                    btInd_sorted = sort(btInd);
                    ci = [btInd_sorted(n*0.025), btInd_sorted(n*0.975)]; % 95% confidence interval
                    MEAN = mean(btInd);
                    bootstrap_summary{1, j} = MEAN;
                    bootstrap_summary{2, j} = ci(1);
                    bootstrap_summary{3, j} = ci(2);
                    if ci(1)<0 && ci(2)>0
                        bootstrap_summary{4, j} = 0;
                    else
                        bootstrap_summary{4, j} = 1;
                    end
                end

                MouseName = se.userData.sessionInfo.MouseName;
                seshDate = se.userData.sessionInfo.seshDate;
                recSite = se.userData.sessionInfo.recSite;

                % Combine with zScore
                new_btInd = [{MouseName} {seshDate} {recSite} {bootstrap_summary} {{drifting_summary}}]; 
                % bootstrap_summary includes mean, lower ci, upper ci, and significant
                % results for each unit
                % drifting_summary: 1 is drifting and 0 is not drifting

                % Concatenation 
                btInds = [btInds; new_btInd];

                clearvars -except seFilePaths btInds MouseName

            end

            % Convert to a table
            VarNames = {'MouseName', 'seshDate', 'recSite', 'btInd', 'drifting'};
            BlockType_Ind = cell2table(btInds,'VariableNames', VarNames);

            % Save btInd summary for each mouse 
            btIndPaths = fullfile('E:\CM_NeuralActivity_Analysis\btInd_atleast4blocks',[MouseName,'_btInd_dfInd']);
            save(btIndPaths, 'BlockType_Ind');
        end
        
        function allConcat(btIndDir)
            % Concatenation 
            btIndPaths = MBrowse.Files(btIndDir);
            btInd_all = table();
            for i = 1: length(btIndPaths)
                load(btIndPaths{i})
                btInd_all = [btInd_all; BlockType_Ind];
            end

            % Save btInx summary for all mice 
            save('E:\CM_NeuralActivity_Analysis\btInd_atleast4blocks\btInd_dfInd', 'btInd_all');
        end
        
        function ranSelect(btInd_tetrode, n)
            % Randomly select n number of units from tetrode dataset 
            % Concatenation
            btInd = [];
            drifting = [];
            for j = 1:size(btInd_tetrode,1) % loop each session
                btInd = [btInd btInd_tetrode{j,4}{1,1}];
                drifting = [drifting btInd_tetrode{j,5}{1,1}]; 
            end
            rng(19); % control random number generation
            randomInd = randperm(size(btInd,2),n); % randomly select n number from 1 tototal number of units
            btInd_selected = btInd(:,randomInd);
            drifting_selected = drifting(:,randomInd);
            recSite = btInd_tetrode.recSite{1};
            btInd_tetrode_selected = [{'Random'} {'Random'} {recSite} {btInd_selected} {{drifting_selected}}]; 
            % Convert to a table
            VarNames = {'MouseName', 'seshDate', 'recSite', 'btInd', 'drifting'};
            BlockType_Ind = cell2table(btInd_tetrode_selected,'VariableNames', VarNames);

            % Save btInd summary for each mouse 
            btIndPaths = fullfile('E:\CM_NeuralActivity_Analysis\btInd_atleast4blocks',...
                ['Random', '_', recSite, '_btInd_dfInd']);
            save(btIndPaths, 'BlockType_Ind');            
        end
        
        function tetrode2siConcat(btIndDir, btInd_SiProbe)
            % Input
            %   btIndDir  where to find block type index of tetrode dataset
            %   btInd_SiProbe block type index of silicone probe dataset
            % Concatenation 
            btIndPaths = MBrowse.Files(btIndDir);
            btInd_all = btInd_SiProbe;
            for i = 1: length(btIndPaths)
                load(btIndPaths{i})
                btInd_all = [btInd_all; BlockType_Ind];
            end

            % Save btInx summary for all mice 
            save('E:\CM_NeuralActivity_Analysis\btInd_atleast4blocks\btInd_dfInd_mixed', 'btInd_all');
        end
  
        function Plot_BlockType(btInd_all, recSite, save_path)
            % Plot btInd distribution for recording sites 
            % session with at least 4 blocks and the last block is longer than 30 trials

            x_all = [];
            y_all = [];
            x_sig = [];
            y_sig = [];
            x_nonsig = [];
            y_nonsig = [];
            perc_sig = [];

            for i = 1: length(recSite)
                t1 = btInd_all(strcmp(btInd_all.recSite, recSite{i}), :); % select recording sites
                for j = 1:size(t1,1) % loop each session
                    t2 = t1{j,4}{1,1}; % get btInd summary
                    t3 = cell2mat(t2(1,:)); % get means
                    df_sig = logical(cell2mat(t1{j,5}{1,1})); % get drifting index
                    new_x = repmat(i,1,length(t3));
                    x_all = [x_all new_x];
                    y_all = [y_all t3];
                    % Get stable units with significant btInd 
                    ind_sig = logical(cell2mat(t2(4, :)));
                    new_y_sig = t3(ind_sig & ~df_sig);
                    new_x_sig = repmat(i,1,length(new_y_sig));
                    x_sig = [x_sig new_x_sig];
                    y_sig = [y_sig new_y_sig];
                    % Get stable units with non significant btInd
                    ind_nonsig = ~ind_sig;
                    new_y_nonsig = t3(ind_nonsig & ~df_sig);
                    new_x_nonsig = repmat(i,1,length(new_y_nonsig));
                    x_nonsig = [x_nonsig new_x_nonsig];
                    y_nonsig = [y_nonsig new_y_nonsig];   
                end
                % percentage of significant units in recording sites
                perc_sig(i) = sum(x_sig(:) == i)/sum(x_all(:) == i);
                xtickNames(i) = strcat(recSite(i),'(', num2str(round(perc_sig(i),2)*100), '%)');
                num_units(i) = sum(x_all(:) == i); % number of units 
            end

            clf
            scatter(x_nonsig,y_nonsig,10,'k', 'jitter','on', 'jitterAmount',0.2)
            hold on
            scatter(x_sig,y_sig,10, 'r', 'filled', 'jitter','on', 'jitterAmount',0.2)
            xlim([0.5 length(recSite)+0.5]);
            ylim([-1.2 1.2]);
            ylabel({'Block type Index'; '(t-v)/(t+v)'})
            xticklabels(xtickNames)
            set(gca,'xtick',1:1:8)
            xtickangle(45)
            
            % Save fig
            btFigPath = fullfile(save_path,'btInd_nondrifting_mixed_somSys.pdf');
            print(btFigPath,'-dpdf','-painters','-loose')
        end
        
        function unitPlot_BlockType(MouseName, seDir)
            
            % Inputs
            %   MouseName    Mouse name
            %   seDir        Directory of session explore
            % Output
            %   Figures of nondrifting units with significant bolck type index
            %   in sessions with at least 4 blocks 
            
            seFilePaths = MBrowse.Files(seDir);
            for j = 1: length(seFilePaths)
                load(seFilePaths{j})

                % Load session information 
                MouseName = se.userData.sessionInfo.MouseName;
                seshDate = se.userData.sessionInfo.seshDate;
                recSite = se.userData.sessionInfo.recSite;

                %%%%%%%%%%%%% Session selection %%%%%%%%%%%%%
                % Selection session with at least 4 blocks
                % Get onset and offset of block trainsition

                % Create index for block numbers 
                trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                             'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
                ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
                        'VariableNames', trialTypes);
                BlockChangeInd = ismember(se.GetColumn('behavValue', 'blockType'), 'Change') ;
                tBlockInd = sum(ttInd{:,1:6},2); 
                vBlockInd = sum(ttInd{:,7:12},2);

                % Get transition offset for each block (automatically reward cue or the first hit, which one comes first)
                % and assign block type to the cahnge trials 
                tBlock_offset = [];
                vBlock_offset = [];
                for t = 2: length(tBlockInd)
                    if BlockChangeInd(t) == 1 && tBlockInd(t-1) ==1
                        tBlockInd(t) = 1;
                        tBlock_offset = [tBlock_offset, t];
                    elseif BlockChangeInd(t) == 1 && vBlockInd(t-1) ==1
                        vBlockInd(t) = 1;
                        vBlock_offset = [vBlock_offset, t];
                    end 
                end

                % Find transition onset (block onset)
                [tBlock_onset,~] = find(diff([0;tBlockInd])==1); % onset = 1; offset = -1
                [vBlock_onset,~] = find(diff([0;vBlockInd])==1);

                % Get how many trials in the last tblock and vBlock
                tBlock_endLength = length(tBlockInd) - tBlock_onset(length(tBlock_onset));
                vBlock_endLength = length(vBlockInd) - vBlock_onset(length(vBlock_onset));

                % Stop if number of blocks are less than 2 
                if length(tBlock_onset) <2 || length(vBlock_onset) <2 
                    continue
                elseif length(tBlock_onset) == 2 && length(vBlock_onset) == 2 
                    % Stop if it is 4 blocks but the last block is shorter than 30 trials 
                    if tBlock_endLength < 30 || vBlock_endLength < 30
                        continue
                    end
                end

                %%%%%%%%%%%%% Get and plot "non-drifting" and "significant block index" units %%%%%%%%%%%%%

                % Assign block numbers within each block type (1st tBlock =1, 2nd tBlock =2...)
                tBlockInd_number = tBlockInd;
                vBlockInd_number = vBlockInd;
                for t = 1: length(tBlockInd_number) % each trial
                    for n = 2: length(tBlock_onset) % each tBlock
                        if t >= tBlock_onset(n) && tBlockInd_number(t) == n-1
                            tBlockInd_number(t) = tBlockInd_number(t) + 1;
                        end
                    end   
                    for m = 2: length(vBlock_onset)
                        if t >= vBlock_onset(m) && vBlockInd_number(t) == m-1
                            vBlockInd_number(t) = vBlockInd_number(t) + 1;
                        end
                    end
                end

                % Get mean baseline activity (firing rate (Hz) for each trials, bin = 1 sec)
                spkRs_mean = CM.Analysis.GetSpikeRate(se, 1, -1, 0); % (Hz)
                baseLine = table2cell(spkRs_mean(:,2:end));
                tBlockInds = logical(repmat(tBlockInd,1,size(baseLine,2)));
                vBlockInds = logical(repmat(vBlockInd,1,size(baseLine,2)));

                tBaseLine = reshape(baseLine(tBlockInds),[sum(tBlockInd), size(baseLine,2)]);
                vBaseLine = reshape(baseLine(vBlockInds),[sum(vBlockInd), size(baseLine,2)]);


                for n = 1: length(tBlock_onset)
                    tBlockInds_block = repmat(ismember(tBlockInd_number,[n]), 1, size(baseLine,2));
                    tBaseLine_block(n) = {reshape(baseLine(tBlockInds_block),[sum(ismember(tBlockInd_number,[n])), size(baseLine,2)])};
                end 
                for m = 1: length(vBlock_onset)
                    vBlockInds_block = repmat(ismember(vBlockInd_number,[m]), 1, size(baseLine,2));
                    vBaseLine_block(m) = {reshape(baseLine(vBlockInds_block),[sum(ismember(vBlockInd_number,[m])), size(baseLine,2)])};
                end 

                % Get mean lick rate (Hz) for each trials, bin = 1 sec
                behavTime_mean = se.ResampleEventTimes('behavTime', -1:1:0);
                rLickRate = behavTime_mean.rLickOnset;
                lLickRate = behavTime_mean.lLickOnset;
                % Get spike data 
                spk = se.GetTable('spikeTime');
                % Get lick data
                rLick = se.GetColumn('behavTime', 'rLickOnset');
                lLick = se.GetColumn('behavTime', 'lLickOnset');

                for i = 1:size(baseLine,2) % unit

                    %%% Select "non-drifting" and "significant block type index" units %%% 
                    % Determinr whether drifting happens (skip if drifting occurs)
                    tBlock_unit = [];
                    vBlock_unit = [];
                    tBlock_grp = [];
                    vBlock_grp = [];

                    % ANOVA: compare baseline activity among different tblocks 
                    for n = 1:length(tBlock_onset)
                        tBlock_unit = [tBlock_unit; tBaseLine_block{1,n}(:,i)]; % Concatenation
                        tBlock_grp = [tBlock_grp; repmat(n,size(tBaseLine_block{1,n},1) ,1)]; % Create group number
                    end
                    % ANOVA: compare baseline activity among different vblocks 
                    for m = 1:length(vBlock_onset)
                        vBlock_unit = [vBlock_unit; vBaseLine_block{1,m}(:,i)]; % Concatenation
                        vBlock_grp = [vBlock_grp; repmat(m,size(vBaseLine_block{1,m},1) ,1)]; % Create group number
                    end

                    p1 = anova1(cell2mat(tBlock_unit), tBlock_grp, 'off');
                    p2 = anova1(cell2mat(vBlock_unit), vBlock_grp, 'off');

                    if p1 < 0.05 || p2 < 0.05
                        continue
                    end

                    % Determine whether baselines differ cross block types
                    n = 10000; % nboot
                    t_bootstat = bootstrp(n,@mean,cell2mat(tBaseLine(:,i)));
                    v_bootstat = bootstrp(n,@mean,cell2mat(vBaseLine(:,i)));
                    btInd = (t_bootstat - v_bootstat)./(t_bootstat + v_bootstat); % Index for block type
                    btInd_sorted = sort(btInd);
                    ci = [btInd_sorted(n*0.025), btInd_sorted(n*0.975)]; % 95% confidence interval
                    MEAN = mean(btInd);

                    % Stop if block type index is not significant
                    if ci(1)<0 && ci(2)>0
                        continue
                    end

                    %%% Plotting %%%
                    s = spk{:,i}; % Select unit

                    % Raster plot (spikes and licks)
                    figure('Position', [0,0, 350, 800]); % set figure size
                    subplot('Position', [0.2 0.1 0.35 0.8])
                    ind = 1:size(s,1);
                    MPlot.PlotRaster(s(ind), ind, 0.5, 'Color', 'k'); hold on
                    MPlot.PlotRaster(rLick(ind), ind, 0.5, 'Color', 'b','LineWidth', 4); hold on
                    MPlot.PlotRaster(lLick(ind), ind, 0.5, 'Color', 'r', 'LineWidth', 4);
                    xlim([-1 0]);
                    ylim([ind(1)-0.5 length(s)+0.5]);
                    xlabel('Time from stimulus onset');
                    ylabel('Trials');
                    % firing rate plot (spikes and licks)
                    subplot('Position', [0.575 0.1 0.25 0.8])
                    smoothed_FR = smoothdata(cell2mat(baseLine(:,i)), 'gaussian', 6);
                    smoothed_rLR = smoothdata(cell2mat(rLickRate), 'gaussian', 6);
                    smoothed_lLR = smoothdata(cell2mat(lLickRate), 'gaussian', 6);
                    y3 = 1:size(s,1);
                    plot(smoothed_FR, y3, 'k', smoothed_rLR, y3, 'b', smoothed_lLR, y3, 'r' ); 
                    xlabel('Firing\lick rate (Hz)');
                    xlim([0 max(smoothed_FR)+1]);
                    ylim([ind(1)-0.5 length(s)+0.5]);
                    set(gca,'yticklabel',{[]});
                    % Block type bar
                    subplot('Position', [0.845 0.1 0.05 0.8])
                    x1 = repmat(0,[size(tBaseLine, 1),1]);
                    y1 = find(tBlockInd == 1);
                    x2 = repmat(0,[size(vBaseLine, 1),1]);
                    y2 = find(vBlockInd == 1);
                    scatter(x1,y1,12,'b','s','filled');
                    hold on
                    scatter(x2,y2,12,'r','s','filled');
                    ylim([ind(1)-0.5 length(s)+0.5]);
                    set(gca,'visible','off')
                    
                    % Silicon probe
                    sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(i),'\_',...
                        num2str(se.userData.sessionInfo.unit_position(i)),'um'],'FontSize',10);
                    % Tetrode
%                     sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(i)],'FontSize',10);     

                    % Save fig
                    figDir = 'E:\CM_NeuralActivity_Analysis\UnitPlot_nondrifting_sigBlockType_licking';
                    cd(figDir);
                    figName = [MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(i)];
                    print(figName,'-dpdf','-painters','-loose');

                end
                clearvars -except seFilePaths
                close all
                
            end  
        end
    end
end