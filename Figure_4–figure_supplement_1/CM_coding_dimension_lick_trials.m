%% Figure_4–figure supplement_1 g-h
% Calculate choice coding dimensions using lick right and lick left trials in respond-to-touch and respond-to-light blocks
% For stimulus coding dimensions, (1) all tactile vs visual trials (including no lick trials) or (2) tactile vs visual trials
% (excluding no lick trials).

load('E:\Cross-Modal_Project\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% setting 
mainDir = 'E:\Cross-Modal_Project\Revision\coding_dimensions';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
B=2;
S=2;
D=3;
bin = 0.01;
rng(15)

timeWindow = [-0.1 0.15]; % normalization window 
analysisOn_Bin = round((timeWindow(1)+1)*(1/bin));
analysisOff_Bin = round((timeWindow(2)+1)*(1/bin)) + 1;
time = timeWindow(1):bin:timeWindow(2)+bin;
time_bw = [time, fliplr(time)];
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

nboot = 1000;
shuffle_conditions = [0 1];

for i = 1: length(recSites)
    u = table2cell(data_all(strcmp(data_all.recSite, recSites{i}), :));
    isSession = cell2mat(cellfun(@(x) size(x,1)>10, u(:,4), 'UniformOutput',false)); % sessions with more than 10 units are included
    u = u(isSession,:);
    for p = 1:nboot
        % resample sessions
        isSession = randi(size(u,1),size(u,1),1);
        resampled_u = u(isSession,:);
        for session=1:size(resampled_u,1)
            Xtrial = resampled_u{session,4}(:,:,:,:,analysisOn_Bin:analysisOff_Bin,:);
            trialNum = resampled_u{session,5};
            % resample neurons and trials
            [resampled_Xtrial, resampled_trialNum] = CM.Statistics.Multi_bootstrp(Xtrial,trialNum,B,S,D);
            
            for q =1:2
                isShuffle = shuffle_conditions(q);
                if isShuffle
                    [X, XtrialShuffle] = CM.Statistics.getShuffledData(resampled_Xtrial,resampled_trialNum);
                else
                    for b = 1:B
                        for s = 1:S
                            for d = 1:D
                                resampled_trialNum_tt = resampled_trialNum(1,b,s,d);
                                resampled_Xtrial(:,b,s,d,:,resampled_trialNum_tt+1:end) = nan;
                            end
                        end
                    end
                    X = nanmean(resampled_Xtrial,6);
                end

                X_mean = nanmean(X(:,:),2);
                X_range = range(X(:,:),2); 
                % normalization
                normalized_X = (X - X_mean)./(X_range + 0.01); 

                % Compute coding dimensions
                for b=1:B
%                     switch b
%                         case 1 % respond-to-touch block
%                             lickDirection = 1; % right lick
%                             block_name = 'respond-to-touch block';
%                         case 2 % respond-to-light block
%                             lickDirection = 2; % left lick
%                             block_name = 'respond-to-light block';
%                     end
                    wt = squeeze(nanmean(normalized_X(:,b,1,1:2,:),4) - nanmean(normalized_X(:,b,2,1:2,:),4)); % tactile stimulus trials - visual stimulus trials
                    CD_stim_mode = mean(wt(:,processingOn_Bin:processingOff_Bin),2);
                    CD_stim_mode = CD_stim_mode/norm(CD_stim_mode);
                    
                    wt = squeeze(nanmean(normalized_X(:,b,1:S,1,:),3) - nanmean(normalized_X(:,b,1:S,2,:),3)); % right lick trials - left lick trials 
                    CD_choice_mode = mean(wt(:,processingOn_Bin:processingOff_Bin),2);
                    CD_choice_mode = CD_choice_mode/norm(CD_choice_mode);

                    CD_recSite(q,session,b,1) = {CD_stim_mode};
                    CD_recSite(q,session,b,2) = {CD_choice_mode};

                    % projections onto coding dimensions
                    projection_on_stimulus_CD(q,session,b,1,1,1:T) = squeeze(normalized_X(:,b,1,1,:))'*CD_stim_mode; % sP_tStim_rLick 
                    projection_on_stimulus_CD(q,session,b,1,2,1:T) = squeeze(normalized_X(:,b,1,2,:))'*CD_stim_mode; % sP_tStim_lLick 
                    projection_on_stimulus_CD(q,session,b,1,3,1:T) = squeeze(normalized_X(:,b,1,3,:))'*CD_stim_mode; % sP_tStim_noLick
                    
                    projection_on_stimulus_CD(q,session,b,2,1,1:T) = squeeze(normalized_X(:,b,2,1,:))'*CD_stim_mode; % sP_vStim_rLick
                    projection_on_stimulus_CD(q,session,b,2,2,1:T) = squeeze(normalized_X(:,b,2,2,:))'*CD_stim_mode; % sP_vStim_lLick
                    projection_on_stimulus_CD(q,session,b,2,3,1:T) = squeeze(normalized_X(:,b,2,3,:))'*CD_stim_mode; % sP_vStim_noLick
                    
                    projection_on_choice_CD(q,session,b,1,1,1:T) = squeeze(normalized_X(:,b,1,1,:))'*CD_choice_mode; % cP_tStim_rLick
                    projection_on_choice_CD(q,session,b,1,2,1:T) = squeeze(normalized_X(:,b,1,2,:))'*CD_choice_mode; % cP_tStim_lLick
                    projection_on_choice_CD(q,session,b,1,3,1:T) = squeeze(normalized_X(:,b,1,3,:))'*CD_choice_mode; % cP_tStim_noLick
                    
                    projection_on_choice_CD(q,session,b,2,1,1:T) = squeeze(normalized_X(:,b,2,1,:))'*CD_choice_mode; % cP_vStim_rLick
                    projection_on_choice_CD(q,session,b,2,2,1:T) = squeeze(normalized_X(:,b,2,2,:))'*CD_choice_mode; % cP_vStim_lLick
                    projection_on_choice_CD(q,session,b,2,3,1:T) = squeeze(normalized_X(:,b,2,3,:))'*CD_choice_mode; % cP_vStim_noLick
                end
                % Absolute value of dot product of stimulus/choice CD in respond-to-touch and respond-to-light blocks
                DP_stim(q,session) = abs(CD_recSite{q,session,1,1}'*CD_recSite{q,session,2,1}); % stimulus dimension
                DP_choice(q,session) = abs(CD_recSite{q,session,1,2}'*CD_recSite{q,session,2,2}); % choice dimension
            end
        end
        % Average across sessions for each sample
        resampled_projection_on_stimulus_CD(1:2,1:B,1:S,1:D,1:T,p) = squeeze(nanmean(projection_on_stimulus_CD,2));
        resampled_projection_on_choice_CD(1:2,1:B,1:S,1:D,1:T,p) = squeeze(nanmean(projection_on_choice_CD,2));
        resampled_DP_stim(1:2,p) = mean(DP_stim,2);
        resampled_DP_choice(1:2,p) = mean(DP_choice,2);
    end
    CD_summary{i,1}= recSites{i};
    CD_summary{i,2}= resampled_projection_on_stimulus_CD;
    CD_summary{i,3}= resampled_projection_on_choice_CD;
    CD_summary{i,4}= resampled_DP_stim;
    CD_summary{i,5}= resampled_DP_choice;
    clearvars resampled_projection_on_stimulus_CD resampled_projection_on_choice_CD resampled_DP_stim resampled_DP_choice
end
cdPaths = fullfile(mainDir,['CD_' num2str(nboot) 'nBoot']);
save(cdPaths, 'CD_summary');
%% Plot
load('E:\Cross-Modal_Project\Revision\coding_dimensions\CD_noLick_trials_for_sCDs\CD_1000nBoot.mat')
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
nboot = 1000;
% Projections onto CD 
colors_SD = {[0 0 1] [0 0.5 1] [0 1 1]; [1 0 0] [1 0 0.5] [1 0 1]};  % rows: tactile stimulus, visual stimulus; columns: right lick, left lick, no lick
block_names = {'Respond-to-touch block','Respond-to-light block'};
CD_names = {'stimulus','choice'};
shuffle_names = {'true' 'shuffled'};
shuffle_colors = {[0 1 1] [0.5 0.5 0.5]};
for q=1:2 % q=1 for true dataset, q=2 for shuffled dataset
    for i = 1: length(recSites)
        figure('Position',[0 0 350 1000])
        for j=1:2 % j=1 for stimulus dimension, j=2 for choice dimension
            resampled_projection = CD_summary{i,j+1};
            for b=1:B
                subplot(4,1,(j-1)*2+b)
                for s=1:S
                    for d=1:D % d=1 for right lick; d=2 for left lick, d=3 for no lick
                        projection_tt = num2cell(squeeze(resampled_projection(q,b,s,d,:,:)),2);
                        projection_sorted = cellfun(@(x) sort(x), projection_tt, 'UniformOutput', false);
                        projection_sorted_noNAN = cellfun(@(x) x(~isnan(x)), projection_sorted, 'UniformOutput', false);
                        bootstrap_summary = cellfun(@(x) [mean(x), x(round(0.025*length(x))), x(round(0.975*length(x)))], projection_sorted_noNAN,'UniformOutput', false);
                        bootstrap_summary = cell2mat(bootstrap_summary); % rows: time bin, columns: mean, lower CI, upper CI  
                        plot(time,bootstrap_summary(:,1),'Color', colors_SD{s,d}); hold on;
                        fill(time_bw, [bootstrap_summary(:,2)', fliplr(bootstrap_summary(:,3)')],...
                            colors_SD{s,d},'FaceAlpha', 0.3, 'Linestyle', 'none');
                    end
                end
                plot(stimWindow, [1.2 1.2],'k','Linewidth',2)
                xlabel('Time (s)')
                xlim(timeWindow);ylim([-0.6 1.2]);
                set(gca, 'box','off','TickDir','out')
                ylabel({block_names{b};'Projection on'; [CD_names{j},' dimension (AU)']})
            end
        end
        sgtitle(recSite_names{i})
        % Save figure  
        FigPath = fullfile(mainDir,['CD_',recSite_names{i},'_bootstrapping', shuffle_names{q}, '.pdf']);
        print(FigPath,'-dpdf','-painters','-loose');
    end
end

% Dot product (true vs shuffled datasets)
figure('Position',[0 0 800 300])
for j=1:2 % j=1 for stimulus dimension, j=2 for choice dimension
    subplot(1,2,j)
    for i = 1: length(recSites)
        for q=1:2 % q=1 for true dataset, q=2 for shuffled dataset
            y= CD_summary{i,j+3}(q,:);
            sorted_y = sort(y);
            mean_y = mean(sorted_y);
            lowCI = sorted_y(round(0.025*nboot));
            highCI = sorted_y(round(0.975*nboot));
            x=i+q*0.2-0.3;
            errorbar(x,round(mean_y,2),round(mean_y-lowCI,2),round(highCI-mean_y,2), 'o','Color', shuffle_colors{q},...
                'MarkerEdgeColor',shuffle_colors{q},'MarkerFaceColor',shuffle_colors{q},...
                'CapSize', 6, 'LineWidth', 0.5); hold on;
            if i==length(recSites) && j == 2
                text(3.5,(0.8-0.05*q),shuffle_names{q},'Color',shuffle_colors{q})
            end
            summary{j,i,q} = [round(mean_y,2), round(lowCI,2), round(highCI,2)];
        end
    end
    ylabel('|f1*f2|')
    xlim([0.5 4.5]);ylim([0 1]);
    yticks(0:0.1:1);yticklabels({0 [] 0.2 [] 0.4 [] 0.6 [] 0.8 [] 1});
    xticks([1:4]);xticklabels(recSite_names);
    set(gca, 'box','off','TickDir','out')
    title(CD_names{j})
end
 % Save figure  
FigPath = fullfile(mainDir,'Dot_product_bootstrapping.pdf');
print(FigPath,'-dpdf','-painters','-loose');   

