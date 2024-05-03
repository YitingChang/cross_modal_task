%% Neural distance
% For each session, trail-averaged tactile hit and tactile correct rejection were used to obtain the top three principal components
% via PCA.
% Activity for each tactile hit and correct rejection trial was projected onto the top-3 PC space. 
% Distance was calculated for each pair of tactile hit and tactile correct rejection at each time point.
% Chang et al. Figure 4b-c

% load data array
load('E:\Data_array\AllTrials_10bin_Smooth\data_array.mat')
% setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\PCA\Distance\Distance_and_pre-stim_state\Trial_activity';
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
B=2;
S=2;
D=3;
bin = 0.01;
timeWindow = [-0.1 0.15];
pre_timeBin = round((timeWindow(1)+1)*(1/bin));
post_timeBin = round((timeWindow(2)+1)*(1/bin)) + 1;
T = post_timeBin - pre_timeBin +1;
initialWindow = [-0.1 0];
initialOn_Bin = (initialWindow(1)-timeWindow(1))*(1/bin) + 1;
initialOff_Bin = (initialWindow(2)-timeWindow(1))*(1/bin);
stimWindow = [0 0.15];
stimOn_Bin = (stimWindow(1)-timeWindow(1))*(1/bin) + 1;
stimOff_Bin = (stimWindow(2)-timeWindow(1))*(1/bin);
componentsToAnalyze = 3;
Colors_tHit = [0 0 1];
Colors_tCR = [0 1 1];

for i = 1: length(recSites)
    u = table2cell(data_all(strcmp(data_all.recSite, recSites{i}), :));
    summary_recSite(1:size(u,1),1:3) = u(1:size(u,1),1:3);
    for session=1:size(u,1)
        Xtrial = u{session,4}(:,:,:,:,pre_timeBin:post_timeBin,:);
        numOfTrials = u{session,5};
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
        Z = (X - X_mean)./ X_range;  
        Ztrial = (Xtrial - X_mean)./ X_range;  
        N = size(Xtrial,1);
        if N > componentsToAnalyze  
            Z_tHit = squeeze(Z(:,1,1,1,:)); 
            Z_tCR = squeeze(Z(:,2,1,3,:)); 
            Z_tHit_tCR = [Z_tHit Z_tCR];
            % singular value decomposition
            [W,~,~] = svd(Z_tHit_tCR, 'econ');
            % reconstruction (trajectories in the PC coordinate system)
            dataDim = size(Z_tHit);
            R_tHit_tCR =  Z_tHit_tCR'*W;
            % explained variance
            explVar.totalVar = sum(sum(Z_tHit_tCR.^2)); % total variance of original data
            explVar.component = sum(R_tHit_tCR.^2)/ explVar.totalVar * 100; % component variance of reconstructed data (%)
            explVar.cumulative = cumsum(sum(R_tHit_tCR.^2)/ explVar.totalVar * 100); % cumulative variance of reconstructed data (%)
            
            numOfTrials_tHit = numOfTrials(1,1,1,1);
            numOfTrials_tCR = numOfTrials(1,2,1,3);
            Ztrial_tHit = squeeze(Ztrial(:,1,1,1,:,:));
            Ztrial_tCR = squeeze(Ztrial(:,2,1,3,:,:));
            for j = 1:numOfTrials_tHit
                Z_E = Ztrial_tHit(:,:,j);
                R_E =  Z_E'*W;
                E2PC = reshape(R_E(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
                E2PC = num2cell(E2PC,1);
                for k=1:numOfTrials_tCR
                    Z_F = Ztrial_tCR(:,:,k);
                    R_F =  Z_F'*W;
                    F2PC = reshape(R_F(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
                    F2PC = num2cell(F2PC,1);
                    % Distance at each time point
                    distance = cellfun(@(x,y) norm(x-y),E2PC,F2PC,'UniformOutput', false);
                    summary_session(j,k) = {cell2mat(distance)};
                    clear distance
                end
            end
            summary_recSite{session,4} = N;   
            summary_recSite{session,5} = summary_session;
            clear summary_session
            for j = 1:numOfTrials_tHit
                Z_E = Ztrial_tHit(:,:,j);
                R_E =  Z_E'*W;
                E2PC = reshape(R_E(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
                E2PC = num2cell(E2PC,1);
                for k=1:numOfTrials_tHit
                    Z_F = Ztrial_tHit(:,:,k);
                    R_F =  Z_F'*W;
                    F2PC = reshape(R_F(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
                    F2PC = num2cell(F2PC,1);
                    % Distance at each time point
                    distance = cellfun(@(x,y) norm(x-y),E2PC,F2PC,'UniformOutput', false);
                    if j==k
                        summary_session{j,k} = repmat(nan,1,T);
                    else
                        summary_session{j,k} = cell2mat(distance);
                    end
                    clear distance
                end
            end
            summary_recSite{session,6} = summary_session;
            clear summary_session
%             % plot
%             R_tHit = Z_tHit'*W;
%             tHit2PC = reshape(R_tHit(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
%             R_tCR = Z_tCR'*W;
%             tCR2PC = reshape(R_tCR(:,1:componentsToAnalyze)', [componentsToAnalyze dataDim(2:end)]);
%             
%             % 3D figure
%             figure('Position',[0 0 400 400])
%             plot3(tHit2PC(1,:),tHit2PC(2,:),tHit2PC(3,:),'Color',Colors_tHit,'LineWidth', 1); hold on;
%             plot3(tCR2PC(1,:),tCR2PC(2,:),tCR2PC(3,:),'Color',Colors_tCR,'LineWidth', 1); hold on;
%             % stimulus onset and offset
%             scatter3(tHit2PC(1,stimOn_Bin),tHit2PC(2,stimOn_Bin),tHit2PC(3,stimOn_Bin),70,'MarkerFaceColor',Colors_tHit,'MarkerEdgeColor','k'); 
%             scatter3(tHit2PC(1,stimOff_Bin),tHit2PC(2,stimOff_Bin),tHit2PC(3,stimOff_Bin),70,'*','MarkerFaceColor',Colors_tHit,'MarkerEdgeColor',Colors_tHit);         
%             scatter3(tCR2PC(1,stimOn_Bin),tCR2PC(2,stimOn_Bin),tCR2PC(3,stimOn_Bin),70,'MarkerFaceColor',Colors_tCR,'MarkerEdgeColor','k'); 
%             scatter3(tCR2PC(1,stimOff_Bin),tCR2PC(2,stimOff_Bin),tCR2PC(3,stimOff_Bin),70,'*','MarkerFaceColor',Colors_tCR,'MarkerEdgeColor',Colors_tCR); 
%             xlabel(['PC1 (' num2str(round(explVar.component(1),1)) '%)'],'FontSize',10);
%             ylabel(['PC2 (' num2str(round(explVar.component(2),1)) '%)'], 'FontSize',10);
%             zlabel(['PC3 (' num2str(round(explVar.component(3),1)) '%)'], 'FontSize',10);
%             title([recSite_names{i} '\_' summary_recSite{session,2}])
%             FigPath = fullfile(mainDir,'Figures_distance_session',['Distance_' recSite_names{i} '_' summary_recSite{session,2} '_v3']);
%             print(FigPath,'-dpdf','-painters','-loose');
            
        else
            summary_recSite{session,4} = N;
            summary_recSite{session,5} = [];
            summary_recSite{session,6} = [];
        end
        clear Xtrial numOfTrials
    end
    summary{i,1} = recSites{i};
    summary{i,2} = summary_recSite;
    clear summary_recSite
end

%% Plot individual brain areas
figure('Position',[0 0 800 200])
x = timeWindow(1):bin:timeWindow(2)+bin;
for i = 1: length(recSites)
    subplot(1,4,i)
    summary_recSite = summary{i,2};
    for session=1:size(summary_recSite,1)
        distance_session = cell2mat(reshape(summary_recSite{session,5},[],1));
        distance_mean = mean(distance_session,1);
        distance_recSite(session,1:T) = distance_mean;
        plot(x,distance_mean,'Color',[0.75 0.75 0.75]); hold on;
    end 
    y = mean(distance_recSite);
    plot(x,y,'k'); hold on;
    plot(stimWindow, [2.5 2.5],'k','Linewidth',1.5)
    xlim(timeWindow)
    xticks(timeWindow(1):0.05:timeWindow(2))
    xticklabels({-0.1 [] 0 [] 0.1 [] 0.2 []})
    ylim([0 2.5])
    xlabel('Time from stimulus onset')
    set(gca, 'box','off','TickDir','out')
    title(recSite_names{i})
    if i==1
        ylabel({'Distance between'; 'tHit and tCR (AU)'})
    end
    clear distance_recSite 
end
% Save figure 
FigPath = fullfile(mainDir,'Distance_tHit_tCR_ind_p0top15');
print(FigPath,'-dpdf','-painters','-loose');

%% Plot by different pre-stimulus states (close vs far)- all brain areas
Colors_recSites = {[0.8 0.95 0.66] [0.6 0.8 0.4] [0.35 0.6 0.2] [0.1 0.1 0.1]};
percentile = 50;
distance_all = [];
ind_all = [];
pre_stim_all = [];
for i = 1: length(recSites)
    summary_recSite = summary{i,2};
    for session=1:size(summary_recSite,1)
        distance_session = cell2mat(reshape(summary_recSite{session,5},[],1));
        distance_mean = mean(distance_session,1);
        pre_stim_distance(session,1) = mean(distance_mean(initialOn_Bin:initialOff_Bin));
        distance_recSite(session,1:T) = distance_mean;
    end 
    
    pre_stim_all = [pre_stim_all; pre_stim_distance];
    distance_all = [distance_all; distance_recSite];
    ind_recSite = repmat(i,size(summary_recSite,1),1);
    ind_all = [ind_all;ind_recSite];
    clear pre_stim_distance distance_recSite 
end

pre_stim_close = pre_stim_all <= prctile(pre_stim_all,percentile);
pre_stim_far = pre_stim_all >= prctile(pre_stim_all,100-percentile);
x = timeWindow(1):bin:timeWindow(2)+bin;

figure('Position',[0 0 600 600])
for session = 1:length(pre_stim_all)
    distance = distance_all(session,:);
    if pre_stim_close(session)
        LineColor = [1 0.75 0.5];
    elseif pre_stim_far(session)
        LineColor = [0.75 1 0.5];
    else
        LineColor = [0.9 0.9 0.9];
    end
%     plot(x,distance, LineStyle,'Color', Colors_recSites{ind_all(session)}); hold on;
    plot(x,distance,'Color',LineColor); hold on;
end
plot(x,mean(distance_all(pre_stim_close,:)),'Color', [0.55 0.1 0.15]);
plot(x,mean(distance_all(pre_stim_far,:)),'Color', [0.1 0.55 0.15]);        


plot(stimWindow, [2.5 2.5],'k','Linewidth',1.5)
xlim(timeWindow)
xticks(timeWindow(1):0.05:timeWindow(2))
xticklabels({-0.1 [] 0 [] 0.1 [] 0.2 []})
ylim([0 2.5])
xlabel('Time from stimulus onset')
ylabel('Distance(AU)')
set(gca, 'box','off','TickDir','out')

% Permutation test
% Determine whether the distance between tHit and tCR trejectories for close and far groups was siginificant, 
% as measured by the Euclidean distance between them
% O'Connor et al., Neuron, 2010 (Fig. S7D)

rng(10)
nShuffle = 1000;
analysisWindow = [0 0.15];
analysisOn_Bin = (analysisWindow(1)-timeWindow(1))*(1/bin) + 1;
analysisOff_Bin = (analysisWindow(2)-timeWindow(1))*(1/bin);

sessionNum_A = sum(pre_stim_close);
sessionNum_B = sum(pre_stim_far);
Xsession_A = distance_all(pre_stim_close,analysisOn_Bin:analysisOff_Bin);
Xsession_B = distance_all(pre_stim_far,analysisOn_Bin:analysisOff_Bin);
X_A = mean(Xsession_A); 
X_B = mean(Xsession_B); 
distance = norm(X_A - X_B);
Xsession_AB = [Xsession_A; Xsession_B];
for shuffle=1:nShuffle
    new_sessionNum = randperm(sessionNum_A+sessionNum_B);
    isA = new_sessionNum <= sessionNum_A;
    Xsession_A_shuffled = Xsession_AB(isA,:);
    Xsession_B_shuffled = Xsession_AB(~isA,:);
    X_A_shuffled = mean(Xsession_A_shuffled); 
    X_B_shuffled = mean(Xsession_B_shuffled); 
    distance_shuffled(shuffle) = norm(X_A_shuffled - X_B_shuffled);
end
pValue = sum(distance_shuffled>=distance)/length(distance_shuffled);% one tailed

title([num2str(percentile) ' vs ' num2str(100-percentile) '\_p = ' num2str(pValue)])


% Save figure 
FigPath = fullfile(mainDir,['Distance_tHit_tCR_' num2str(percentile) ' vs ' num2str(100-percentile) '_p0top15']);
print(FigPath,'-dpdf','-painters','-loose');



