% Coding dimension
% Find coding dimensions for stimulus and choice during respond-to-touch and respond-to-light rules
% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\Coding_dimensions';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
B=2;
S=2;
D=3;
bin = 0.01;

timeWindow = [-0.1 0.15]; % normalization window 
analysisOn_Bin = round((timeWindow(1)+1)*(1/bin));
analysisOff_Bin = round((timeWindow(2)+1)*(1/bin)) + 1;
time = timeWindow(1):bin:timeWindow(2)+bin;
T = length(time);

initialWindow = [timeWindow(1) 0]; % subspace overlap window (before stimulus onset)
initialOn_Bin = (initialWindow(1)-timeWindow(1))*(1/bin) + 1;
initialOff_Bin = (initialWindow(2)-timeWindow(1))*(1/bin);

processingWindow = [0 timeWindow(2)];  % subspace overlap window (after stimulus onset)
processingOn_Bin = (processingWindow(1)-timeWindow(1))*(1/bin) + 1;
processingOff_Bin = (processingWindow(2)-timeWindow(1))*(1/bin);

stimWindow = [0 0.15];
stimOn_Bin = (stimWindow(1)-timeWindow(1))*(1/bin) + 1;
stimOff_Bin = (stimWindow(2)-timeWindow(1))*(1/bin);

isShuffle = 1;
for i = 1: length(recSites)
    u = table2cell(data_all(strcmp(data_all.recSite, recSites{i}), :));
    isSession = cell2mat(cellfun(@(x) size(x,1)>10, u(:,4), 'UniformOutput',false)); % sessions with more than 10 units are included
    u = u(isSession,:);
    for session=1:size(u,1)
        seshDate = u{session,2};
        Xtrial = u{session,4}(:,:,:,:,analysisOn_Bin:analysisOff_Bin,:);
        numOfTrials = u{session,5};
        if isShuffle
            [X, XtrialShuffle] = CM.Statistics.getShuffledData(Xtrial,numOfTrials);
        else
            for b = 1:B
                for s = 1:S
                    for d = 1:D
                        numOfTrials_tt = numOfTrials(1,b,s,d);
                        Xtrial(:,b,s,d,:,numOfTrials_tt+1:end) = nan;
                    end
                end
            end
            X = nanmean(Xtrial,6);
        end
        X_mean = nanmean(X(:,:),2);
        X_range = range(X(:,:),2); 
        % normalization
        normalized_X = (X - X_mean)./(X_range + 0.01); 
        
        % Compute coding dimensions
%         figure('Position',[0 0 600 600])
        for b=1:B
            switch b
                case 1 % respond-to-touch block
                    lickDirection = 1; % right lick
                    block_name = 'respond-to-touch block';
                case 2 % respond-to-light block
                    lickDirection = 2; % left lick
                    block_name = 'respond-to-light block';
            end
            wt = squeeze((normalized_X(:,b,1,lickDirection,:)+ normalized_X(:,b,1,3,:))/2 ...
                - (normalized_X(:,b,2,lickDirection,:)+ normalized_X(:,b,2,3,:))/2);
            CD_stim_mode = mean(wt(:,processingOn_Bin:processingOff_Bin),2);
            CD_stim_mode = CD_stim_mode/norm(CD_stim_mode);

            wt = squeeze((normalized_X(:,b,1,lickDirection,:)+ normalized_X(:,b,2,lickDirection,:))/2 ...
                - (normalized_X(:,b,1,3,:)+ normalized_X(:,b,2,3,:))/2);
            CD_choice_mode = mean(wt(:,processingOn_Bin:processingOff_Bin),2);
            CD_choice_mode = CD_choice_mode/norm(CD_choice_mode);
            
            CD_recSite(session,b,1) = {CD_stim_mode};
            CD_recSite(session,b,2) = {CD_choice_mode};
            
            % projections onto coding dimensions
            projection_on_stimulus_CD(session,b,1,1,1:T) = squeeze(normalized_X(:,b,1,lickDirection,:))'*CD_stim_mode; % sP_tLick 
            projection_on_stimulus_CD(session,b,1,2,1:T) = squeeze(normalized_X(:,b,1,3,:))'*CD_stim_mode; % sP_tNoLick
            projection_on_stimulus_CD(session,b,2,1,1:T) = squeeze(normalized_X(:,b,2,lickDirection,:))'*CD_stim_mode; % sP_vLick
            projection_on_stimulus_CD(session,b,2,2,1:T) = squeeze(normalized_X(:,b,2,3,:))'*CD_stim_mode; % sP_vNoLick
            projection_on_choice_CD(session,b,1,1,1:T) = squeeze(normalized_X(:,b,1,lickDirection,:))'*CD_choice_mode; % cP_tLick
            projection_on_choice_CD(session,b,1,2,1:T) = squeeze(normalized_X(:,b,1,3,:))'*CD_choice_mode; % cP_tNoLick
            projection_on_choice_CD(session,b,2,1,1:T) = squeeze(normalized_X(:,b,2,lickDirection,:))'*CD_choice_mode; % cP_vLick
            projection_on_choice_CD(session,b,2,2,1:T) = squeeze(normalized_X(:,b,2,3,:))'*CD_choice_mode; % cP_vNoLick
            
%             subplot(2,2,(b-1)*2+1)
%             plot(time,sP_tLick,'-b',time,sP_tNoLick,'--b',time,sP_vLick,'-r',time,sP_vNoLick,'--r')
%             set(gca, 'box','off','TickDir','out')
%             xlabel('Time (s)')
%             ylabel({block_name; 'Projection onto stimulus dimension'})
% 
%             subplot(2,2,(b-1)*2+2)
%             plot(time,cP_tLick,'-b',time,cP_tNoLick,'--b',time,cP_vLick,'-r',time,cP_vNoLick,'--r')
%             set(gca, 'box','off','TickDir','out')
%             xlabel('Time (s)')
%             ylabel({block_name; 'Projection onto choice dimension'})

        end
%         sgtitle([recSite_names{i}, '\_', seshDate])
%         % Save figure 
%         FigPath = fullfile(mainDir,[recSite_names{i}, '_', seshDate, '.pdf']);
%         print(FigPath,'-dpdf','-painters','-loose');
        % Absolute value of dot product of stimulus/choice CD in respond-to-touch and respond-to-light blocks
        DP_stim_recSite(session) = abs(CD_recSite{session,1,1}'*CD_recSite{session,2,1}); % stimulus dimension
        DP_choice_recSite(session) = abs(CD_recSite{session,1,2}'*CD_recSite{session,2,2}); % choice dimension
    end
    CD_summary{i,1}= recSites{i};
    CD_summary{i,2}= CD_recSite;
    CD_summary{i,3}= DP_stim_recSite;
    CD_summary{i,4}= DP_choice_recSite;
    CD_summary{i,5}= projection_on_stimulus_CD;
    CD_summary{i,6}= projection_on_choice_CD;
    clearvars CD_recSite DP_stim_recSite DP_choice_recSite projection_on_stimulus_CD projection_on_choice_CD
end

%% Plot 
% Plot projection
colors_SB = {[0 0 1] [0 1 1]; [1 0 0] [1 0 1]};  % rows: tactile stimulus, visual stimulus; columns: touch block, light block 
lines_D = {'-' '--'};
for i = 1: length(recSites)
    figure('Position',[0 0 300 1000])
    for j=1:2 % j=1 for stimulus dimension, j=2 for choice dimension
        projection_session = CD_summary{i,j+4};
        for b=1:B
            subplot(4,1,(j-1)*2+b)
            for s=1:S
                for d=1:2 % lick and no lick- right lick in tBlock, left in vBlock
                    y = squeeze(mean(projection_session(:,b,s,d,:),1));
                    plot(time,y,'Color',colors_SB{s,b},'LineStyle',lines_D{d}); hold on;
                end
            end
            plot(stimWindow, [1.5 1.5],'k','Linewidth',2)
            xlabel('Time (s)')
            xlim(timeWindow);ylim([-0.5 1.5]);
            set(gca, 'box','off','TickDir','out')
            switch j
                case 1
                    ylabel({'Projection on'; 'stimulus dimension (AU)'})
                case 2
                    ylabel({'Projection on'; 'choice dimension (AU)'})
            end
        end
    end
    sgtitle(recSite_names{i})
    % Save figure  
    if isShuffle
        FigPath = fullfile(mainDir,['CD_',recSite_names{i},'_shuffled.pdf']);
    else
        FigPath = fullfile(mainDir,['CD_',recSite_names{i},'.pdf']);
    end
    print(FigPath,'-dpdf','-painters','-loose');
end


% Plot dot product
figure('Position',[0 0 800 300])
for j=1:2 % j=1 for stimulus dimension, j=2 for choice dimension
    subplot(1,2,j)
    for i = 1: length(recSites)
        y= CD_summary{i,j+2};
        mean_y = mean(y);
        std_y = std(y);
        sem_y = std_y/length(y)^0.5;
        x = repmat(i-0.2,length(y),1);
        scatter(x,y,'b');hold on;
        errorbar(i+0.2,mean_y,sem_y, 'o', 'Color', 'k', 'MarkerSize', 8,'MarkerEdgeColor','k',...
            'MarkerFaceColor', 'b','CapSize', 8, 'LineWidth', 1); 
    end
    ylabel('|f1*f2|')
    xlim([0.5 4.5]);ylim([0 1]);
    yticks(0:0.1:1);yticklabels({0 [] 0.2 [] 0.4 [] 0.6 [] 0.8 [] 1});
    xticks([1:4]);xticklabels(recSite_names);
    set(gca, 'box','off','TickDir','out')
    switch j
        case 1
            title('stimulus')
        case 2
            title('choice')
    end
end
 % Save figure  
if isShuffle
    FigPath = fullfile(mainDir,['Dot_product_of_coding_dimensions_shuffled.pdf']);
else
    FigPath = fullfile(mainDir,['Dot_product_of_coding_dimensions.pdf']);
end
print(FigPath,'-dpdf','-painters','-loose');     
        
% Projections for all sessions
colors_SD = {[0 0 1] [0 1 1]; [1 0 0] [1 0 1]};  % rows: tactile stimulus, visual stimulus; columns: lick, no lick
block_names = {'Respond-to-touch block','Respond-to-light block'};
CD_names = {'stimulus','choice'};
for i = 1: length(recSites)
    figure('Position',[0 0 350 1000])
    for j=1:2 % j=1 for stimulus dimension, j=2 for choice dimension
        projection_session = CD_summary{i,j+4};
        for b=1:B
            subplot(4,1,(j-1)*2+b)
            for s=1:S
                for d=1:2 % lick and no lick- right lick in tBlock, left in vBlock
                    y = squeeze(projection_session(:,b,s,d,:));
                    x = repmat(time,size(y,1),1);
                    plot(x',y','Color',colors_SD{s,d}); hold on;
                end
            end
            plot(stimWindow, [2 2],'k','Linewidth',2)
            xlabel('Time (s)')
            xlim(timeWindow);ylim([-1 2]);
            set(gca, 'box','off','TickDir','out')
            ylabel({block_names{b};'Projection on'; [CD_names{j},' dimension (AU)']})
        end
    end
    sgtitle(recSite_names{i})
    % Save figure  
    FigPath = fullfile(mainDir,['CD_',recSite_names{i},'_sessions.pdf']);
    print(FigPath,'-dpdf','-painters','-loose');
end

