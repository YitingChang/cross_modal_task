%% Plot correct trial types grouped by stim types
% Plot nondrifting unit (at least 4 blocks and the last block is longer than 30 trials) 
MouseName = 'YT084';
MainDir = 'F:\';
seDir = [MainDir, MouseName, '\MSessionExplorer_', MouseName];
seFilePaths = MBrowse.Files(seDir);
 

for j = 1: length(seFilePaths)
    load(seFilePaths{j})
    
    %%%%%%%% Create index for block numbers %%%%%%%%
    trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
                 'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
    ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
            'VariableNames', trialTypes);

    BlockChangeInd = ismember(se.GetColumn('behavValue', 'blockType'), 'Change') ;

    tBlockInd = sum(ttInd{:,1:6},2); 
    vBlockInd = sum(ttInd{:,7:12},2);
    % Assign block type to the cahnge trials  
    for t = 2: length(tBlockInd)
        if BlockChangeInd(t) == 1 && tBlockInd(t-1) ==1
            tBlockInd(t) = 1;
        elseif BlockChangeInd(t) == 1 && vBlockInd(t-1) ==1
            vBlockInd(t) = 1;
        end 
    end

    % Find block onset 
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
    
    %%%%%%%% Get spike rates of trial types for each unit
    bin = 0.025;
    spkRs = CM.Analysis.GetSpikeRate(se, bin, -1, 1);
    spkRs_unit = spkRs(:,2:end);%remove time column
    spkRate_tt = CM.Analysis.GetUnitSummary(se, spkRs_unit);
    
    
    %%%%%%%% Get drifting and block index for each unit %%%%%%%%
    spkRs_base = CM.Analysis.GetSpikeRate(se, 1, -1, 0); % (Hz)
    baseLine = table2cell(spkRs_base(:,2:end));
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
       
    
    for i = 1:size(spkRate_tt.TTT{1,1},2) % unit
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
        if ci(1)<0 && ci(2)>0
            btInd_sig = 0;
        else
            btInd_sig = 1;
        end
        
        % Get spike times
        spk = se.GetTable('spikeTime');
%         trialTypes = {'TTT', 'TTV', 'TTN', 'TVT', 'TVV', 'TVN',...
%                           'VTT', 'VTV', 'VTN', 'VVT', 'VVV', 'VVN'};
%         ttInd = array2table(se.GetColumn('behavValue', trialTypes),...
%             'VariableNames', trialTypes);
        s = spk{:,i}; % Select unit
        % tactile trials
        s1 = s(ttInd{:,1}); % Select trial types: TTT
%         s2 = s(ttInd{:,3}); % TTN
%         s3 = s(ttInd{:,8}); % VTV
        s4 = s(ttInd{:,9}); % VTN
        
        % visual trials
        s5 = s(ttInd{:,11}); % VVV
%         s6 = s(ttInd{:,12}); % VVN
%         s7 = s(ttInd{:,4}); % TVT
        s8 = s(ttInd{:,6}); % TVN
        
        % reaction time
        lickOn = se.GetColumn('behavTime', 'firstLick');
        ttt_lickOn = num2cell(lickOn(ttInd{:,1})); % TTT
%         vtv_lickOn = num2cell(lickOn(ttInd{:,8})); % VTV
        vvv_lickOn = num2cell(lickOn(ttInd{:,11})); % VVV 
%         tvt_lickOn = num2cell(lickOn(ttInd{:,4})); % TVT
        s1 = [s1 ttt_lickOn];
%         s3 = [s3 vtv_lickOn];
        s5 = [s5 vvv_lickOn];
%         s7 = [s7 tvt_lickOn];
        
        sorted_s1 = sortrows(s1,2);
%         sorted_s3 = sortrows(s3,2);
        sorted_s5 = sortrows(s5,2);
%         sorted_s7 = sortrows(s7,2);

%         s_tStim = [sorted_s1(:,1); s2; sorted_s3(:,1); s4];
%         s_vStim = [sorted_s5(:,1); s6; sorted_s7(:,1); s8];
        s_tStim = [sorted_s1(:,1); s4];
        s_vStim = [sorted_s5(:,1); s8];
        
        % Load session information 
        MouseName = se.userData.sessionInfo.MouseName;
        seshDate = se.userData.sessionInfo.seshDate;
        recSite = se.userData.sessionInfo.recSite;
        
        % trial number
        tt_trialNum = size(sorted_s1,1);
        vt_trialNum = size(s4,1);
        tStim_trialNum = tt_trialNum + vt_trialNum;
        
        vv_trialNum = size(sorted_s5,1);
        tv_trialNum = size(s8,1);
        vStim_trialNum = vv_trialNum + tv_trialNum;
        
%         tt_trialNum = size(sorted_s1,1)+size(s2,1);
%         vt_trialNum = size(sorted_s3,1)+size(s4,1);
%         tStim_trialNum = tt_trialNum + vt_trialNum;
%         
%         vv_trialNum = size(sorted_s5,1)+size(s6,1);
%         tv_trialNum = size(sorted_s7,1)+size(s8,1);
%         vStim_trialNum = vv_trialNum + tv_trialNum;
        
        % block length
        [M,I] = max([tStim_trialNum vStim_trialNum]);
        if I == 1
            tStim_barRatio = tStim_trialNum/tStim_trialNum;
            vStim_barRatio = vStim_trialNum/tStim_trialNum;
        else
            tStim_barRatio = tStim_trialNum/vStim_trialNum;
            vStim_barRatio = vStim_trialNum/vStim_trialNum;
        end
        
        % Plot
        % tStim raster plot
        figure('Position', [0,0, 480, 600]);clf % set figure size
        subplot('Position', [0.2 0.35 0.25 0.5]) 
        ind = 1:size(s_tStim,1);
        MPlot.PlotRaster(s_tStim(ind), ind, 0.5, 'Color', 'k');
        ttt_trailNum = (1:size(sorted_s1,1)).';
        scatter(cell2mat(sorted_s1(:,2)), ttt_trailNum, 8, 'filled', 'b');
%         vtv_trailNum = (tt_trialNum+1:tt_trialNum+size(sorted_s3,1)).';
%         scatter(cell2mat(sorted_s3(:,2)), vtv_trailNum, 8, 'filled', 'g');
        xlim([-0.25 1]);
        ylim([0.5 M+0.5]);
        set(gca,'xticklabel',{[]})
        ylabel('Trials');
        title('Tactile trials', 'FontSize',8);
        box off
        
        % Block type bar
        subplot('Position', [0.453 0.35 0.05 0.5*tStim_barRatio])
        x1 = repmat(0,[tt_trialNum,1]);
        y1 = (1:tt_trialNum).';
        x2 = repmat(0,[vt_trialNum,1]);
        y2 = (tt_trialNum+1:tt_trialNum+vt_trialNum).';
%         scatter(x1,y1,12,[1 0.68 0.4],'s','filled');
        scatter(x1,y1,12,[0.3 0.7 0],'s','filled');

        hold on
%         scatter(x2,y2,12,[1 0.4 0.16],'s','filled');
        scatter(x2,y2,12,[0.5 0.9 0.5],'s','filled');

        ylim([0.5 tt_trialNum+vt_trialNum+0.5]);
        set(gca,'visible','off')
        
        % vStim raster plot
        subplot('Position', [0.6 0.35 0.25 0.5])
        ind = 1:size(s_vStim,1);
        MPlot.PlotRaster(s_vStim(ind), ind, 0.5, 'Color', 'k');
        vvv_trailNum = (1:size(sorted_s5,1)).';
        scatter(cell2mat(sorted_s5(:,2)), vvv_trailNum, 8,'filled', 'r'); 
%         tvt_trailNum = (vv_trialNum+1:vv_trialNum+size(sorted_s7,1)).';
%         scatter(cell2mat(sorted_s7(:,2)), tvt_trailNum, 8, 'filled', 'g');
        xlim([-0.25 1]);
        ylim([ind(1)-0.5 M+0.5]);
        set(gca,'xticklabel',{[]})
        ylabel('Trials');
        title('Visual trials', 'FontSize',8);
        box off
        
        % Block type bar
        subplot('Position', [0.853 0.35 0.05 0.5*vStim_barRatio])
        
        x1 = repmat(0,[vv_trialNum,1]);
        y1 = (1:vv_trialNum).';
        x2 = repmat(0,[tv_trialNum,1]);
        y2 = (vv_trialNum+1:vv_trialNum+tv_trialNum).';
%         scatter(x1,y1,12,[1 0.4 0.16],'s','filled');
        scatter(x1,y1,12,[0.5 0.9 0.5],'s','filled');

        hold on
        scatter(x2,y2,12,[0.3 0.7 0],'s','filled');
%         scatter(x2,y2,12,[1 0.68 0.4],'s','filled');
        ylim([0.5 vv_trialNum+tv_trialNum+0.5]);
        set(gca,'visible','off')
        
        % Average trace
        time = spkRs{1,1}{1,1}.';
        subplot('Position', [0.2 0.1 0.25 0.2])
        plot(time, spkRate_tt.TTT{1,1}(:,i), 'b',time, spkRate_tt.VTN{1,1}(:,i),'--b');
        hold on
        CI_TTT = [spkRate_tt.TTT{3,1}(:,i).', fliplr(spkRate_tt.TTT{4,1}(:,i).')];
%         CI_TTN = [spkRate_tt.TTN{3,1}(:,i).', fliplr(spkRate_tt.TTN{4,1}(:,i).')];
%         CI_VTV = [spkRate_tt.VTV{3,1}(:,i).', fliplr(spkRate_tt.VTV{4,1}(:,i).')];
        CI_VTN = [spkRate_tt.VTN{3,1}(:,i).', fliplr(spkRate_tt.VTN{4,1}(:,i).')];
        time_bw = [time, fliplr(time)];
        fill(time_bw, CI_TTT, 'b','FaceAlpha', 0.5, 'Linestyle', 'none');
        fill(time_bw, CI_VTN, 'b', 'FaceAlpha', 0.3, 'Linestyle', 'none');

        xlim([-0.25 1]);
        ylim([-0.5 4]);
        xlabel('Time from stim onset(s)');
        ylabel('Spike rate (/25ms)');
        hold off
        
        subplot('Position', [0.6 0.1 0.25 0.2])
        plot(time, spkRate_tt.VVV{1,1}(:,i),'r',time, spkRate_tt.TVN{1,1}(:,i),'--r');
        hold on
        CI_VVV = [spkRate_tt.VVV{3,1}(:,i).', fliplr(spkRate_tt.VVV{4,1}(:,i).')];
%         CI_VVN = [spkRate_tt.VVN{3,1}(:,i).', fliplr(spkRate_tt.VVN{4,1}(:,i).')];
%         CI_TVT = [spkRate_tt.TVT{3,1}(:,i).', fliplr(spkRate_tt.TVT{4,1}(:,i).')];
        CI_TVN = [spkRate_tt.TVN{3,1}(:,i).', fliplr(spkRate_tt.TVN{4,1}(:,i).')];
        time_bw = [time, fliplr(time)];
        fill(time_bw, CI_VVV, 'r', 'FaceAlpha', 0.5, 'Linestyle', 'none');
        fill(time_bw, CI_TVN, 'r','FaceAlpha', 0.3, 'Linestyle', 'none');

        xlim([-0.25 1]);
        ylim([-0.5 4]);
        xlabel('Time from stim onset(s)');
        ylabel('Spike rate (/25ms)');
        hold off

        sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(i),...
        '\_btSig\_', num2str(btInd_sig)],'FontSize',10);
        
%         sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(i),'\_',...
%             num2str(se.userData.sessionInfo.unit_position(i)),'um\_btSig\_', num2str(btInd_sig)],'FontSize',10);
        % Save fig
        figDir = 'E:\CM_NeuralActivity_Analysis\UnitPlot_btInd_nondrifting_correctTrials';
        cd(figDir);
        figName = [MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(i)];
        print(figName,'-dpdf','-painters','-loose');

        %saveFigurePDF(fig,...
        %    ['E:\UnitPlot\',MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(i), '.pdf']);
        
    end
    clearvars -except seFilePaths
    close all
end