%% Unit Plot
% Plot nondrifting unit (at least 4 blocks and the last block is longer than 30 trials) 
MouseName = 'EF0114';
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
        s1 = s(ttInd{:,1}); % Select trial types: TTT
        s2 = s(ttInd{:,6}); % TVN
        s3 = s(ttInd{:,11}); % VVV
        s4 = s(ttInd{:,9}); % VTN

        % reaction time
        lickOn = se.GetColumn('behavTime', 'firstLick');
        ttt_lickOn = num2cell(lickOn(ttInd{:,1})); % TTT
        vvv_lickOn = num2cell(lickOn(ttInd{:,11})); % VVV 
        s1 = [s1 ttt_lickOn];
        s3 = [s3 vvv_lickOn];

        sorted_s1 = sortrows(s1,2);
        sorted_s3 = sortrows(s3,2);

        s_tBlock = [sorted_s1(:,1); s2];
        s_vBlock = [sorted_s3(:,1); s4];

        % Load session information 
        MouseName = se.userData.sessionInfo.MouseName;
        seshDate = se.userData.sessionInfo.seshDate;
        recSite = se.userData.sessionInfo.recSite;

        % Plot
        % tBlock raster plot
        figure('Position', [0,0, 480, 600]);clf % set figure size
        subplot('Position', [0.2 0.35 0.25 0.5])
        ind = 1:size(s_tBlock,1);
        MPlot.PlotRaster(s_tBlock(ind), ind, 0.5, 'Color', 'k');
        t_trailNum = (1:length(sorted_s1)).';
        scatter(cell2mat(sorted_s1(:,2)), t_trailNum, 8, 'filled', 'b');  
        xlim([-1 1]);
        ylim([ind(1)-0.5 max([length(s_tBlock) length(s_vBlock)])+0.5]);
        set(gca,'xticklabel',{[]})
        ylabel('Trials');
        hold on
        rectangle('Position',[1 0 0.1 length(s1)], 'EdgeColor', 'b')
        rectangle('Position',[1 length(s1) 0.1 length(s2)], 'EdgeColor', 'r')
        title('Tactile Block', 'FontSize',8);
        box off

        
        % vBlock raster plot
        subplot('Position', [0.6 0.35 0.25 0.5])
        ind = 1:size(s_vBlock,1);
        MPlot.PlotRaster(s_vBlock(ind), ind, 0.5, 'Color', 'k');
        v_trailNum = (1:length(sorted_s3)).';
        scatter(cell2mat(sorted_s3(:,2)), v_trailNum, 8,'filled', 'r');  
        xlim([-1 1]);
        ylim([ind(1)-0.5 max([length(s_tBlock) length(s_vBlock)])+0.5]);
        set(gca,'xticklabel',{[]})
        ylabel('Trials');
        hold on
        rectangle('Position',[1 0 0.1 length(s3)], 'EdgeColor', 'r')
        rectangle('Position',[1 length(s3) 0.1 length(s4)], 'EdgeColor', 'b')
        title('Visual Block', 'FontSize',8);
        box off

        % Average trace
        time = spkRs{1,1}{1,1}.';
        subplot('Position', [0.2 0.1 0.25 0.2])
        plot(time, spkRate_tt.TTT{1,1}(:,i), 'b',time, spkRate_tt.TVN{1,1}(:,i),'r:');
        hold on
        CI_TTT = [spkRate_tt.TTT{3,1}(:,i).', fliplr(spkRate_tt.TTT{4,1}(:,i).')];
        CI_TVN = [spkRate_tt.TVN{3,1}(:,i).', fliplr(spkRate_tt.TVN{4,1}(:,i).')];
        time_bw = [time, fliplr(time)];
        fill(time_bw, CI_TTT, 'b', time_bw, CI_TVN, 'r','FaceAlpha', 0.3, 'Linestyle', 'none');
        xlim([-1 1]);
        ylim([-0.5 4]);
        xlabel('Time from stim onset(s)');
        ylabel('Spike rate (/25ms)');
        hold off

        subplot('Position', [0.6 0.1 0.25 0.2])
        plot(time, spkRate_tt.VVV{1,1}(:,i),'r',time,spkRate_tt.VTN{1,1}(:,i),'b:');
        hold on
        CI_VVV = [spkRate_tt.VVV{3,1}(:,i).', fliplr(spkRate_tt.VVV{4,1}(:,i).')];
        CI_VTN = [spkRate_tt.VTN{3,1}(:,i).', fliplr(spkRate_tt.VTN{4,1}(:,i).')];
        time_bw = [time, fliplr(time)];
        fill(time_bw, CI_VVV, 'r', time_bw, CI_VTN, 'b','FaceAlpha', 0.3, 'Linestyle', 'none');
        xlim([-1 1]);
        ylim([-0.5 4]);
        xlabel('Time from stim onset(s)');
        ylabel('Spike rate (/25ms)');
        hold off

        sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(i),...
        '\_btSig\_', num2str(btInd_sig)],'FontSize',10);
        
%         sgtitle([MouseName, '\_', seshDate, '\_', recSite,'\_unit\_',num2str(i),'\_',...
%             num2str(se.userData.sessionInfo.unit_position(i)),'um\_btSig\_', num2str(btInd_sig)],'FontSize',10);
        % Save fig
        figDir = 'E:\CM_NeuralActivity_Analysis\UnitPlot_btInd_nondrifting';
        cd(figDir);
        figName = [MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(i)];
        print(figName,'-dpdf','-painters','-loose');

        %saveFigurePDF(fig,...
        %    ['E:\UnitPlot\',MouseName,'_', seshDate, '_', recSite, '_unit_', num2str(i), '.pdf']);
        
    end
    clearvars -except seFilePaths
    close all
end

