%% PCA 
%load zScore data
save_path = 'E:\CM_NeuralActivity_Analysis\PCA\PCA_figures_20Cosyne';
% Select trial types and time window for pca 
correctTrials = ['TTT', 'VTN', 'VVV', 'TVN'];
tt_numb = [4, 12, 14, 9];
tt_colors = ['b','c', 'r', 'm']; 

tw = [-0.1 0.5];
% zScore data window: -1 ~ 2.5 sec, bin = 25 ms
pre_timeBin = (tw(1)+1)*1000/25 + 1;
post_timeBin = (tw(2)+1)*1000/25;


% Select zScore based on brain regions of interest
recSite = {'left S1', 'left S2', 'left wM2','left AMM','left ALM', 'right ALM', 'right V1', 'left DLS'};
for i = 1: length(recSite)
    pcaData = [];
    u = zScore_all(strcmp(zScore_all.recSite, recSite{i}), :);

    for tt = 1:length(tt_numb)
        v = [];
        for ss = 1:height(u)
            new_v = u{ss, tt_numb(tt)}{1,1}(pre_timeBin:post_timeBin,:); 
            v = [v new_v];
        end
        pcaData = [pcaData; v];
    end
    
    pcaData = pcaData.';

    [coeff,score,latent,tsquared,explained,mu] = pca(pcaData);
    PC1 = coeff(:,1);
    PC2 = coeff(:,2);

    % plotting
    figure('Position', [0,0, 800, 800]); % set figure size
    len = post_timeBin - pre_timeBin +1;
    stimOn_Bin = (0-tw(1))*1000/25; 
    stimOff_Bin = (0.15-tw(1))*1000/25;
    PC1_StimOn = [PC1(stimOn_Bin), PC1(stimOn_Bin+1*len), PC1(stimOn_Bin+2*len), PC1(stimOn_Bin+3*len)];
    PC2_StimOn = [PC2(stimOn_Bin), PC2(stimOn_Bin+1*len), PC2(stimOn_Bin+2*len), PC2(stimOn_Bin+3*len)];
    PC1_StimOff = [PC1(stimOff_Bin), PC1(stimOff_Bin+1*len), PC1(stimOff_Bin+2*len), PC1(stimOff_Bin+3*len)];
    PC2_StimOff = [PC2(stimOff_Bin), PC2(stimOff_Bin+1*len), PC2(stimOff_Bin+2*len), PC2(stimOff_Bin+3*len)];
    
    hold on
    for tt_plot = 1:length(tt_numb)
        plot(PC1((tt_plot-1)*len+1:(tt_plot)*len), PC2((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot),'linewidth', 1)
        plot(PC1((tt_plot-1)*len+stimOn_Bin:(tt_plot-1)*len+stimOff_Bin),...
            PC2((tt_plot-1)*len+stimOn_Bin:(tt_plot-1)*len+stimOff_Bin), tt_colors(tt_plot),...
            'linewidth', 2.5);
        scatter(PC1_StimOn(tt_plot), PC2_StimOn(tt_plot),60,tt_colors(tt_plot),'filled', 'MarkerEdgeColor', 'k');
        %scatter(PC1_StimOff(tt_plot), PC2_StimOff(tt_plot),40,'s',tt_colors(tt_plot), 'filled');
    end
     
    hold off

    xlabel(['PC1(',num2str(round(explained(1),1)),'%)'],'FontSize',12);
    ylabel(['PC2(',num2str(round(explained(2),1)),'%)'],'FontSize',12);
    title([recSite{i}], 'FontSize',16);
    
    fig = gcf;
    PCAFigPath = fullfile(save_path,[recSite{i},'_CorrectTrials_PCA.pdf']);
    print(PCAFigPath,'-dpdf','-painters','-loose');
end
    
%% Plot PC space
%load zScore data
save_path = 'E:\CM_NeuralActivity_Analysis\PCA\PCA_figures_20Cosyne';
% Select trial types and time window for pca 
correctTrials = ['TTT', 'VTN', 'VVV', 'TVN'];
tt_numb = [4, 12, 14, 9];
tt_colors = ['b','c', 'r', 'm']; 

tw = [-0.1 0.5];
% zScore time window: -0.1~0.5 sec, bin = 25 ms
pre_timeBin = (tw(1)+1)*1000/25 + 1;
post_timeBin = (tw(2)+1)*1000/25;


% Select zScore based on brain regions of interest
recSite = {'left S1', 'left S2', 'left wM2','left AMM','left ALM', 'right ALM', 'right V1', 'left DLS'};
for i = 1: length(recSite)
    pcaData = [];
    u = zScore_all(strcmp(zScore_all.recSite, recSite{i}), :);

    for tt = 1:length(tt_numb)
        v = [];
        for ss = 1:height(u)
            new_v = u{ss, tt_numb(tt)}{1,1}(pre_timeBin:post_timeBin,:); 
            v = [v new_v];
        end
        pcaData = [pcaData; v];
    end
    
    pcaData = pcaData.';

    [coeff,score,latent,tsquared,explained,mu] = pca(pcaData);
    PC1 = coeff(:,1);
    PC2 = coeff(:,2);

    % plotting
    figure('Position', [0,0, 800, 600]); % set figure size
    len = post_timeBin - pre_timeBin +1;
    stimOn_Bin = (0-tw(1))*1000/25; 
    stimOff_Bin = (0.15-tw(1))*1000/25;
    
    stimOn_Bins = stimOn_Bin:len:len*length(tt_numb);
    stimOff_Bins = stimOff_Bin:len:len*length(tt_numb);
    for tt_plot = 1:length(tt_numb)
        subplot(2,1,1)
        x = ((tt_plot -1)*len +1):1:((tt_plot)*len);
        plot(x,PC1((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot))
        xline(stimOn_Bins(tt_plot), '--');
        xline(stimOff_Bins(tt_plot), '--');
        hold on
        ylabel(['PC1(',num2str(round(explained(1),1)),'%)'],'FontSize',12);
        xlabel('Time');
        set(gca,'xticklabel',{[]})
    
        subplot(2,1,2)
        plot(x,PC2((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot))
        xline(stimOn_Bins(tt_plot), '--'); 
        xline(stimOff_Bins(tt_plot), '--');
        hold on
        ylabel(['PC2(',num2str(round(explained(2),1)),'%)'],'FontSize',12);
        xlabel('Time');
        set(gca,'xticklabel',{[]})
    end  
    sgtitle(recSite(i))
    
    fig = gcf;
    PCAFigPath = fullfile(save_path,[recSite{i},'_CorrectTrials_PC_space.pdf']);
    print(PCAFigPath,'-dpdf','-painters','-loose');
end
    



%% PCA 
% pooled all recording sites
recSite = {'left S1', 'left S2', 'left wM2','left AMM','left ALM', 'right ALM'};
u = zScore_all(ismember(zScore_all.recSite, recSite), :);

pcaData = [];

for tt = 1:length(tt_numb)
    v = [];
    for ss = 1:height(u)
        new_v = u{ss, tt_numb(tt)}{1,1}(pre_timeBin:post_timeBin,:); 
        v = [v new_v];
    end
    pcaData = [pcaData; v];
end

pcaData = pcaData.';

[coeff,score,latent,tsquared,explained,mu] = pca(pcaData);
PC1 = coeff(:,1);
PC2 = coeff(:,2);

% plotting
figure('Position', [0,0, 800, 800]); % set figure size
len = post_timeBin - pre_timeBin +1;
stimOn_Bin = (0-tw(1))*1000/25; 
stimOff_Bin = (0.15-tw(1))*1000/25;
PC1_StimOn = [PC1(stimOn_Bin), PC1(stimOn_Bin+1*len), PC1(stimOn_Bin+2*len), PC1(stimOn_Bin+3*len)];
PC2_StimOn = [PC2(stimOn_Bin), PC2(stimOn_Bin+1*len), PC2(stimOn_Bin+2*len), PC2(stimOn_Bin+3*len)];
PC1_StimOff = [PC1(stimOff_Bin), PC1(stimOff_Bin+1*len), PC1(stimOff_Bin+2*len), PC1(stimOff_Bin+3*len)];
PC2_StimOff = [PC2(stimOff_Bin), PC2(stimOff_Bin+1*len), PC2(stimOff_Bin+2*len), PC2(stimOff_Bin+3*len)];

hold on
for tt_plot = 1:length(tt_numb)
    plot(PC1((tt_plot-1)*len+1:(tt_plot)*len), PC2((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot),'linewidth', 1)
    plot(PC1((tt_plot-1)*len+stimOn_Bin:(tt_plot-1)*len+stimOff_Bin),...
        PC2((tt_plot-1)*len+stimOn_Bin:(tt_plot-1)*len+stimOff_Bin), tt_colors(tt_plot),...
        'linewidth', 2.5);
    scatter(PC1_StimOn(tt_plot), PC2_StimOn(tt_plot),60,tt_colors(tt_plot),'filled', 'MarkerEdgeColor', 'k');
    %scatter(PC1_StimOff(tt_plot), PC2_StimOff(tt_plot),40,'s',tt_colors(tt_plot), 'filled');
end

hold off

xlabel(['PC1(',num2str(round(explained(1),1)),'%)'],'FontSize',12);
ylabel(['PC2(',num2str(round(explained(2),1)),'%)'],'FontSize',12);
title([recSite{i}], 'FontSize',16);

fig = gcf;
PCAFigPath = fullfile(save_path,'allRecSite_CorrectTrials_PCA.pdf');
print(PCAFigPath,'-dpdf','-painters','-loose');
% How to diassociate recSite later?

%% PCA for baseline activity 
%load zScore data
save_path = 'E:\CM_NeuralActivity_Analysis\PCA\PCA_figures_20Cosyne';
% Select trial types and time window for pca 
correctTrials = ['TTT', 'VTN', 'VVV', 'TVN'];
tt_numb = [4, 12, 14, 9];
tt_colors = ['b','c', 'r', 'm']; 

tw = [-1 0];
% zScore data window: -1 ~ 2.5 sec, bin = 25 ms
pre_timeBin = (tw(1)+1)*1000/25 + 1;
post_timeBin = (tw(2)+1)*1000/25;


% Select zScore based on brain regions of interest
recSite = {'left S1', 'left S2', 'left wM2','left AMM','left ALM', 'right ALM', 'right V1', 'left DLS'};
for i = 1: length(recSite)
    pcaData = [];
    u = zScore_all(strcmp(zScore_all.recSite, recSite{i}), :);

    for tt = 1:length(tt_numb)
        v = [];
        for ss = 1:height(u)
            new_v = u{ss, tt_numb(tt)}{1,1}(pre_timeBin:post_timeBin,:); 
            v = [v new_v];
        end
        pcaData = [pcaData; v];
    end
    
    pcaData = pcaData.';

    [coeff,score,latent,tsquared,explained,mu] = pca(pcaData);
    PC1 = coeff(:,1);
    PC2 = coeff(:,2);

    % plotting
    figure('Position', [0,0, 800, 800]); % set figure size
    len = post_timeBin - pre_timeBin +1;
    
    hold on
    for tt_plot = 1:length(tt_numb)
        plot(PC1((tt_plot-1)*len+1:(tt_plot)*len), PC2((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot),'linewidth', 1)
    end
    hold off

    xlabel(['PC1(',num2str(round(explained(1),1)),'%)'],'FontSize',12);
    ylabel(['PC2(',num2str(round(explained(2),1)),'%)'],'FontSize',12);
    title([recSite{i}], 'FontSize',16);
    
    fig = gcf;
    PCAFigPath = fullfile(save_path,[recSite{i},'_CorrectTrials_PCA_baseline.pdf']);
    print(PCAFigPath,'-dpdf','-painters','-loose');
end

%% Plot PC space for baseline
%load zScore data
save_path = 'E:\CM_NeuralActivity_Analysis\PCA\PCA_figures_20Cosyne';
% Select trial types and time window for pca 
correctTrials = ['TTT', 'VTN', 'VVV', 'TVN'];
tt_numb = [4, 12, 14, 9];
tt_colors = ['b','c', 'r', 'm']; 

tw = [-1 0];
% zScore data window: -1 ~ 2.5 sec, bin = 25 ms
pre_timeBin = (tw(1)+1)*1000/25 + 1;
post_timeBin = (tw(2)+1)*1000/25;

% Select zScore based on brain regions of interest
recSite = {'left S1', 'left S2', 'left wM2','left AMM','left ALM', 'right ALM', 'right V1', 'left DLS'};
for i = 1: length(recSite)
    pcaData = [];
    u = zScore_all(strcmp(zScore_all.recSite, recSite{i}), :);

    for tt = 1:length(tt_numb)
        v = [];
        for ss = 1:height(u)
            new_v = u{ss, tt_numb(tt)}{1,1}(pre_timeBin:post_timeBin,:); 
            v = [v new_v];
        end
        pcaData = [pcaData; v];
    end
    
    pcaData = pcaData.';

    [coeff,score,latent,tsquared,explained,mu] = pca(pcaData);
    PC1 = coeff(:,1);
    PC2 = coeff(:,2);

    % plotting
    figure('Position', [0,0, 800, 600]); % set figure size
    len = post_timeBin - pre_timeBin +1;
    
    for tt_plot = 1:length(tt_numb)
        subplot(2,1,1)
        x = ((tt_plot -1)*len +1):1:((tt_plot)*len);
        plot(x,PC1((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot))
        
        hold on
        ylabel(['PC1(',num2str(round(explained(1),1)),'%)'],'FontSize',12);
        xlabel('Time');
        set(gca,'xticklabel',{[]})
    
        subplot(2,1,2)
        plot(x,PC2((tt_plot-1)*len+1:(tt_plot)*len),tt_colors(tt_plot))
        hold on
        ylabel(['PC2(',num2str(round(explained(2),1)),'%)'],'FontSize',12);
        xlabel('Time');
        set(gca,'xticklabel',{[]})
    end  
    sgtitle(recSite(i))
    
    fig = gcf;
    PCAFigPath = fullfile(save_path,[recSite{i},'_CorrectTrials_PC_space_baseline.pdf']);
    print(PCAFigPath,'-dpdf','-painters','-loose');
end
    