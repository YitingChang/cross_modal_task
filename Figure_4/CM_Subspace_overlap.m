%% Subspace overlap for VTN and control TTT 

% For each session, subspace overlap is calculated by average tHit and tCR. 
% no session concatenation
% Dimentionality reduction analysis is not suitable for a sessions with few units. 
% To control this, only sessions with more than 10 units are included and the number of principal components for analysis is 3.
% Normalization (response - mean)/range

% Touch-evoked activity in tHit trials was used to obtain principal components. 
% Pre-stimulus and touch-evoked activity in tCR trials was then projected in that PC space. 

% Chang et al., Figure 4 

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\PCA\Subspace_overlap\subspace_overlap_control_TTT_stim_window';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
B=2;
S=2;
D=3;
bin = 0.01;

timeWindow = [-0.1 0.15]; % normalization window 
analysisOn_Bin = round((timeWindow(1)+1)*(1/bin));
analysisOff_Bin = round((timeWindow(2)+1)*(1/bin)) + 1;

initialWindow = [timeWindow(1) 0]; % subspace overlap window (before stimulus onset)
initialOn_Bin = (initialWindow(1)-timeWindow(1))*(1/bin) + 1;
initialOff_Bin = (initialWindow(2)-timeWindow(1))*(1/bin);

processingWindow = [0 timeWindow(2)];  % subspace overlap window (after stimulus onset)
processingOn_Bin = (processingWindow(1)-timeWindow(1))*(1/bin) + 1;
processingOff_Bin = (processingWindow(2)-timeWindow(1))*(1/bin);

stimWindow = [0 0.15];
stimOn_Bin = (stimWindow(1)-timeWindow(1))*(1/bin) + 1;
stimOff_Bin = (stimWindow(2)-timeWindow(1))*(1/bin);

componentsToAnalyze = 3;

% conversion map
tt_Es = {[1 1 1]}; % reference: refTTT (tHit)
tt_Fs = {[1 1 1] [2 1 3]}; % comparison: controlTTT VTN (tCR)
% tt_Fs = {[1 1 1] [2 1 2]}; % comparison: controlTTT VTV
Colors_E = {[0 0 1]};
Colors_F2E = {[0.5 0.5 0.5] [0.5 0.2 0.55]};
LineStyles_E = {'-'};
LineStyles_F2E = {'-'};
conType = {'TTT2TTT' 'VTN2TTT'};
nCV = 2;

for i = 1: length(recSites)
    u = table2cell(data_all(strcmp(data_all.recSite, recSites{i}), :));
    isSession = cell2mat(cellfun(@(x) size(x,1)>10, u(:,4), 'UniformOutput',false)); % sessions with more than 10 units are included
    u = u(isSession,:);
    for session=1:size(u,1)
        seshDate = u{session,2};
        Xtrial = u{session,4}(:,:,:,:,analysisOn_Bin:analysisOff_Bin,:);
        numOfTrials = u{session,5};
        N = size(Xtrial,1);
        for b = 1:B
            for s = 1:S
                for d = 1:D
                    numOfTrials_tt = numOfTrials(1,b,s,d);
                    Xtrial(:,b,s,d,:,numOfTrials_tt+1:end) = nan;
                end
            end
        end

        X = nanmean(Xtrial,6);
        X_mean = nanmean(X(:,:),2);
        X_range = range(X(:,:),2); 

        for ref = 1:length(tt_Es)
            tt_E = tt_Es{ref};
            % The post-stimulus-onset window for calculating PCs
            Xtrial2 = Xtrial(:,:,:,:,processingOn_Bin:processingOff_Bin,:);
            T =  size(Xtrial2,5); 
            %%%%%%%%% refernce and comaprison groups %%%%%%%%%
            for n=1:N
                numOfTrials_E = numOfTrials(n,tt_E(1),tt_E(2),tt_E(3));
                Xtrial_n_E = squeeze(Xtrial2(n,tt_E(1),tt_E(2),tt_E(3),:,1:numOfTrials_E));
                isComp = randperm(numOfTrials_E) <= round(numOfTrials_E/nCV,0);
                X_E_comp(n,1:T) = mean(Xtrial_n_E(:,isComp),2);
                X_E_ref(n,1:T) = mean(Xtrial_n_E(:,~isComp),2);
            end

            %%%%%%%%% PCA on the reference group of trial type E %%%%%%%%%
            % normalization
            Z_E_ref = (X_E_ref- X_mean)./(X_range + 5);   

            % singular value decomposition
            [W_E,~,~] = svd(Z_E_ref, 'econ'); 

            for pp = 1:2 % pre-stimulus-onset and post-stimulus-onset
                if pp == 1 % Calculate subspace overlap during pre-stimulus-onset window
                    X2 = X(:,:,:,:,initialOn_Bin:initialOff_Bin);
                elseif pp == 2 % Calculate subspace overlap during post-stimulus-onset window
                    X2 = X(:,:,:,:,processingOn_Bin:processingOff_Bin);
                end
                for comp = 1:length(tt_Fs)
                    tt_F = tt_Fs{comp};
                    %%%%%%%%% PCA on F %%%%%%%%%
                    if pp==2 && sum(tt_F == tt_E)==3 % the comparison group of trial type E during the post-stimulus-onset window
                        Z_F = (X_E_comp- X_mean)./(X_range + 0.01);  
                    else % trial type F and trial type E during the pre-stimulus-onset window
                        X_F = squeeze(X2(:,tt_F(1),tt_F(2),tt_F(3),:));
                        Z_F = (X_F - X_mean)./(X_range + 0.01); 
                    end
                    % dataDim = size(Z_F);
                    % singular value decomposition
                    [W_F,~,~] = svd(Z_F, 'econ'); 
                    % reconstruction onto the top N PC space (trajectories in the PC coordinate system)
                    W_F_topPC = W_F(:,1:componentsToAnalyze);
                    R_F =  Z_F'*W_F_topPC;

                    %%%%%%%%% Map F onto E PC space %%%%%%%%%
                    % reconstruction onto the top N PC space
                    W_E_topPC = W_E(:,1:componentsToAnalyze);
                    R_F2E =  Z_F'*W_E_topPC;
                    % Z_F2E = reshape(R_F2E(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
                    
                    % subspace overlap
                    subspace_overlap = cumsum(sum(R_F2E.^2))./cumsum(sum(R_F.^2));
                    subspace_overlap_recSite(session, pp, comp) = subspace_overlap(componentsToAnalyze);
                end
                clearvars X2 
            end 
            clearvars X_E_comp X_E_ref
        end
        clearvars Xtrial numOfTrials
    end
    subspace_overlap_summary{i,1}= recSites{i};
    subspace_overlap_summary{i,2}= subspace_overlap_recSite;
    clearvars subspace_overlap_recSite
end

%% Plot subspace overlaps for VTN and control TTT
% figure('Position', [0 0 600 800]); 
% for pp = 1:2 % pre-stimulus-onset vs post-stimulus-onset 
%     subplot(2,1,pp)
%     for i = 1: length(recSites)
%         for comp = 1:length(tt_Fs)
%             y= subspace_overlap_summary{i,2}(:,pp,comp);
%             mean_y = mean(y);
%             std_y = std(y);
%             sem_y = std_y/length(y)^0.5;
%             x = repmat(i+(comp-1)*0.3-0.15,length(y),1);
%             scatter(x,y,'MarkerEdgeColor',Colors_F2E{comp});hold on;
% %             errorbar(i+(comp-1)*0.3,mean_y,sem_y, 'o', 'Color', 'k', 'MarkerSize', 8,'MarkerEdgeColor','k',...
% %                 'MarkerFaceColor', Colors_F2E{comp},'CapSize', 8, 'LineWidth', 1); 
%         end
%         NumOfSessions = length(y);
%         x_session = [i-0.15, i+0.15];
%         for s=1:NumOfSessions
%             y_session = squeeze(subspace_overlap_summary{i,2}(s,pp,:));
%             plot(x_session,y_session,'Color',[0.7 0.7 0.7])
%         end     
%     end
%     xlim([0.5 4.5]);ylim([0 1]);
%     yticks(0:0.1:1);yticklabels({0 [] 0.2 [] 0.4 [] 0.6 [] 0.8 [] 1});
%     xticks([1:4]);xticklabels(recSite_names);
%     if pp == 1 % Calculate subspace overlap during pre-stimulus-onset window
%         timeWindow = initialWindow;
%     elseif pp == 2 % Calculate subspace overlap during post-stimulus-onset window
%         timeWindow = processingWindow;
%     end
%             
%     ylabel({'Subspace overlap';[ '(' num2str(timeWindow(1)) ' to ' num2str(timeWindow(2)) ' s)']})
%     set(gca, 'box','off','TickDir','out')
% end

%  % Save figure  
% FigPath = fullfile(mainDir,['Subspace_overlap_control_TTT_post_VTN.pdf']);
% print(FigPath,'-dpdf','-painters','-loose');

% plot pre- and post-stimulus onset subspace overlap in one figure 
figure('Position', [0 0 600 400]); 
for i = 1: length(recSites)
    % post-stimulus onset subspace overlap (TTT and VTN)
    for comp = 1:length(tt_Fs)
        y= subspace_overlap_summary{i,2}(:,2,comp);
%         mean_y = mean(y);
%         std_y = std(y);
%         sem_y = std_y/length(y)^0.5;
        x = repmat(i+(comp-1)*0.6-0.3,length(y),1);
        scatter(x,y,'MarkerEdgeColor',Colors_F2E{comp});hold on;
%             errorbar(i+(comp-1)*0.3,mean_y,sem_y, 'o', 'Color', 'k', 'MarkerSize', 8,'MarkerEdgeColor','k',...
%                 'MarkerFaceColor', Colors_F2E{comp},'CapSize', 8, 'LineWidth', 1); 
    end
    % pre-stimulus onset subspace overlap (VTN)
    y= subspace_overlap_summary{i,2}(:,1,2);
%     mean_y = mean(y);
%     std_y = std(y);
%     sem_y = std_y/length(y)^0.5;
    x = repmat(i,length(y),1);
    scatter(x,y,'MarkerEdgeColor',[0.8 0.5 0.85]);hold on;

    NumOfSessions = length(y);
    x_session = [i-0.3, i, i+0.3];
    for s=1:NumOfSessions
        y_session = [subspace_overlap_summary{i,2}(s,2,1),subspace_overlap_summary{i,2}(s,1,2),...
            subspace_overlap_summary{i,2}(s,2,2)]; % post-stim+TTT, pre-stim+VTN, post-stim+VTN
        plot(x_session,y_session,'Color',[0.7 0.7 0.7])
    end     
end
xlim([0.5 4.5]);ylim([0 1]);
yticks(0:0.1:1);yticklabels({0 [] 0.2 [] 0.4 [] 0.6 [] 0.8 [] 1});
xticks([1:4]);xticklabels(recSite_names);          
ylabel('Subspace overlap')
set(gca, 'box','off','TickDir','out')
 % Save figure  
FigPath = fullfile(mainDir,['Subspace_overlap_control_TTT_pre_post_VTN.pdf']);
print(FigPath,'-dpdf','-painters','-loose');


% pre/post VTN vs control TTT
for i = 1: length(recSites) 
    tHit = subspace_overlap_summary{i,2}(:,2,1); % control TTT (post-stimulus-onset)
    for pp = 1:2 % pre-stimulus-onset vs post-stimulus-onset 
        tCR = subspace_overlap_summary{i,2}(:,pp,2);
        [h,p,ci,stats] = ttest(tHit,tCR); % paired-t test
        so_tHit_tCR{i,pp,1}= recSites{i};
        so_tHit_tCR{i,pp,2}= mean(tHit)- mean(tCR);
        so_tHit_tCR{i,pp,3}= p;
    end
end

% post VTN (comparison between brain areas)
% for i = 1: length(recSites) 
%     tCR_A = subspace_overlap_summary{i,2}(:,2,2); 
%     for j = (i+1):length(recSites) 
%         tCR_B = subspace_overlap_summary{j,2}(:,2,2); 
%         [h,p,ci,stats] = ttest2(tCR_A,tCR_B); % two-sample t test
%         so_tCR{i,j,1}= recSites{i};
%         so_tCR{i,j,2}= recSites{j};
%         so_tCR{i,j,3}= mean(tCR_A)- mean(tCR_B);
%         so_tCR{i,j,4}= p;
%     end
% end

%% Plot subspace overlap for pre- vs post-stimulus-onset VTN
% individual brain areas
Colors_recSites = {[0.8 0.95 0.66] [0.6 0.8 0.4] [0.35 0.6 0.2] [0.1 0.1 0.1]};
figure('Position',[0 0 800 800])
for i = 1: length(recSites)
    x = subspace_overlap_summary{i,2}(:,1,2); % before stimulus onset
    y = subspace_overlap_summary{i,2}(:,2,2);% after stimulus onset
    subplot(2,2,i)
    scatter(x,y,40,'filled','MarkerEdgeColor',Colors_recSites{i},'MarkerFaceColor',Colors_recSites{i});
    xlim([0 1])
    ylim([0 1])
    xlabel({'Subspace overlap before stimulus onset';[ '(' num2str(initialWindow(1)) ' to ' num2str(initialWindow(2)) ' s)']})
    ylabel({'Subspace overlap after stimulus onset';[ '(' num2str(processingWindow(1)) ' to ' num2str(processingWindow(2)) ' s)']})
    title(recSites{i})
    lsline
end
% Save figure 
FigPath = fullfile(mainDir,'subspace_overlap_pre_vs_post_VTN_ind.pdf');
print(FigPath,'-dpdf','-painters','-loose');


% all brain areas
figure('Position',[0 0 400 400])
X = [];
Y = [];
for i = 1: length(recSites)
    x = subspace_overlap_summary{i,2}(:,1,2); % before stimulus onset
    y = subspace_overlap_summary{i,2}(:,2,2);% after stimulus onset
    X = [X;x];
    Y=[Y;y];
end   
scatter(X,Y,40,'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth', 0.1);hold on;
ls = lsline;hold on;
ls.Color = 'k';
[F,P] = corrcoef(X,Y);
title(['r = ' num2str(round(F(1,2),2)) '; p = ' num2str(round(P(1,2),6))])

% lable recSites
for i = 1: length(recSites)
    x = subspace_overlap_summary{i,2}(:,1,2); % before stimulus onset
    y = subspace_overlap_summary{i,2}(:,2,2);% after stimulus onset
    scatter(x,y,40,'filled','MarkerEdgeColor',Colors_recSites{i},'MarkerFaceColor',Colors_recSites{i});hold on;
end 

set(gca, 'box','off','TickDir','out')
xlabel({'Subspace overlap before stimulus onset';[ '(' num2str(initialWindow(1)) ' to ' num2str(initialWindow(2)) ' s)']})
ylabel({'Subspace overlap after stimulus onset';[ '(' num2str(processingWindow(1)) ' to ' num2str(processingWindow(2)) ' s)']})
xticks(0:0.2:1)
yticks(0:0.2:1)
xlim([0 1])
ylim([0 1])

% Save figure 
FigPath = fullfile(mainDir,'subspace_overlap_pre_vs_post_VTN_all.pdf');
print(FigPath,'-dpdf','-painters','-loose');

%% Plot subspace overlaps for VTN based on their pre-stimulus distance
% % pre-stimulus distance results: summary
% % Get close and far groups
% percentile = 50;
% pre_stim_all = [];
% for i = 1: length(recSites)
%     summary_recSite = summary{i,2};
%     for session=1:size(summary_recSite,1)
%         distance_session = cell2mat(reshape(summary_recSite{session,5},[],1));
%         distance_mean = mean(distance_session,1);
%         pre_stim_distance(session,1) = mean(distance_mean(initialOn_Bin:initialOff_Bin));
%     end 
%     
%     pre_stim_all = [pre_stim_all; pre_stim_distance];
%     clear pre_stim_distance 
% end
% 
% pre_stim_close = pre_stim_all <= prctile(pre_stim_all,percentile);
% pre_stim_far = pre_stim_all >= prctile(pre_stim_all,100-percentile);
% 
% session_start = 1;
% figure('Position', [0 0 600 400]); 
% so_close = [];
% so_far = [];
% for i = 1: length(recSites)
%     distance_summary_recSite = summary{i,2};
%     isSession = cell2mat(cellfun(@(x) x>10, distance_summary_recSite(:,4), 'UniformOutput',false));
%     session_end = session_start + length(isSession) -1;
%     pre_stim_recSite = pre_stim_all(session_start:session_end);
%     pre_stim_close_recSite = pre_stim_close(session_start:session_end);
%     session_start = session_start + sum(isSession);
%     is_close = pre_stim_close_recSite(isSession);
%     
%     y= subspace_overlap_summary{i,2}(:,2,2);
%     x = repmat(i,length(y),1);
%     scatter(x,y,'MarkerEdgeColor','b');hold on;
%     scatter(x(is_close),y(is_close),'MarkerEdgeColor','r');
%     
%     so_close = [so_close; y(is_close)];
%     so_far = [so_far; y(~is_close)];
% end
% 
% [h,p,ci,stats] = ttest2(so_close, so_far);
% 
% xlim([0.5 4.5]);ylim([0 1]);
% yticks(0:0.1:1);yticklabels({0 [] 0.2 [] 0.4 [] 0.6 [] 0.8 [] 1});
% xticks([1:4]);xticklabels(recSite_names);
% ylabel({'Subspace overlap';[ '(' num2str(processingWindow(1)) ' to ' num2str(processingWindow(2)) ' s)']})
% set(gca, 'box','off','TickDir','out')
% 
% % Save figure  
% FigPath = fullfile(mainDir,['Subspace_overlap_post_VTN_close_far_groups.pdf']);
% print(FigPath,'-dpdf','-painters','-loose');
% 
% 
% session_start = 1;
% X=[];
% Y=[];
% figure('Position', [0 0 600 400]); 
% for i = 1: length(recSites)
%     distance_summary_recSite = summary{i,2};
%     isSession = cell2mat(cellfun(@(x) x>10, distance_summary_recSite(:,4), 'UniformOutput',false));
%     session_end = session_start + length(isSession) -1;
%     pre_stim_recSite = pre_stim_all(session_start:session_end);
%     pre_stim_close_recSite = pre_stim_close(session_start:session_end);
%     session_start = session_start + sum(isSession);
%     is_close = pre_stim_close_recSite(isSession);
%     
%     y= subspace_overlap_summary{i,2}(:,2,2);
%     x = pre_stim_recSite(isSession);
%     X = [X;x];
%     Y=[Y;y];
% end
% scatter(X,Y,40,'MarkerEdgeColor',[0.5 0.5 0.5], 'LineWidth', 0.1);hold on;
% ls = lsline;hold on;
% ls.Color = 'k';
% [F,P] = corrcoef(X,Y);
% title(['r = ' num2str(round(F(1,2),2)) '; p = ' num2str(round(P(1,2),6))])
% ylim([0 1]);
% yticks(0:0.1:1);yticklabels({0 [] 0.2 [] 0.4 [] 0.6 [] 0.8 [] 1});
% xlabel('Pre-stimulus distance')
% ylabel('Subspace overlap')
% set(gca, 'box','off','TickDir','out')
% 
% 
% % Save figure  
% FigPath = fullfile(mainDir,['Subspace_overlap_vs_preStim_distance.pdf']);
% print(FigPath,'-dpdf','-painters','-loose');
