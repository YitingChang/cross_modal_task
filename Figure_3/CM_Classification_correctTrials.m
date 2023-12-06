%% Classification of initial states of correct trials 
% Trial-by-trial analysis
% For each session, all correct trials were used to calculate classification accuracy. 
% 
% load data array
load('E:\Data_array\AllTrials_100bin_noSmooth\data_array.mat')
% Import setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\CorrectTrials\Trial_by_trial';
B=2;
S=2;
D=3;
bin = 0.1;
timeWindow = [-1 2.5];
analysisWindow = [-1 0.3];
binStart =round((analysisWindow(1)-timeWindow(1))/bin+1);
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);
% initialBin = (0+1)*(1/bin); % -0.1-0s
% initialBin = (0.1+1)*(1/bin); % 0-0.1s
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
tts = {[1 1 1] [1 2 3] [2 1 3] [2 2 2]}; % trial types for classification: TTT (tHit), TVN (vCR), VTN (tCR), VVV (vHit) 
component_names = {'Block' 'Stimulus' 'Decision'};
shuffle_names = {'true' 'shuffled'};
shuffle_colors = {[0 1 1] [0.5 0.5 0.5]};
nCV = 10;
nBoot = 100;
rng(22)
NumTrees = 500;

%% True and shuffled classification accuracy (session)
% Built-in cross-validation. 
% Chang et al., Fig 3e and S3
for binNum= binStart:binEnd
    for i = 1: length(recSites)
        u = table2cell(data_all(strcmp(data_all.recSite, recSites{i}), :));
        for session=1:size(u,1)
%             Xtrial = u{session,4}(:,:,:,:,initialBin,:);
            Xtrial = u{session,4}(:,:,:,:,binNum,:);
            numOfTrials = u{session,5};
            N = size(Xtrial,1); % number of neurons
            %%%%%%%%% Normalization %%%%%%%%
            Xtrial_2D = [];
            for b = 1:B
                for s = 1:S
                    for d = 1:D
                        numOfTrials_tt = numOfTrials(1,b,s,d);
                        Xtrial_tt = squeeze(Xtrial(:,b,s,d,:,1:numOfTrials_tt)); 
                        Xtrial_2D = [Xtrial_2D Xtrial_tt];
                    end
                end
            end
            X_mean = mean(Xtrial_2D,2);
            X_std = std(Xtrial_2D,0,2);
            Ztrial = (Xtrial - X_mean)./(X_std + 0.01); % adding 0.01 to std to prevent the denominator is zero  

            %%%%%%%%% Cross Validation %%%%%%%%
            for shuffle=1:2 % 1: no shuffling; 2: shuffling
                for comp = 1:length(component_names)
                    %%%%%%%%% Concatenation %%%%%%%%
                    isClass = [];
                    Ztrial_2D = [];
                    for tt = 1:length(tts)
                        numOfTrials_tt = numOfTrials(1,tts{tt}(1),tts{tt}(2),tts{tt}(3)); 
                        Ztrial_tt = squeeze(Ztrial(:,tts{tt}(1),tts{tt}(2),tts{tt}(3),1:numOfTrials_tt));
                        Ztrial_2D = [Ztrial_2D Ztrial_tt];
                        if comp ==3 
                            if tts{tt}(comp) == 2 % merge left and right lick
                                isClass_tt = ones(numOfTrials_tt,1);
                            else
                                isClass_tt = repmat(tts{tt}(comp),[numOfTrials_tt 1]);
                            end
                        else
                            isClass_tt = repmat(tts{tt}(comp),[numOfTrials_tt 1]);
                        end
                        isClass = [isClass; isClass_tt];
                        clear isClass_tt Ztrial_tt
                    end 
                    newOrder = randperm(length(isClass));
                    isClass = isClass(newOrder);
                    if shuffle == 1 % true
                        Ztrial_2D = Ztrial_2D(:,newOrder);
                    else % shuffle
                        newOrder2 = randperm(length(isClass));
                        Ztrial_2D = Ztrial_2D(:,newOrder2);
                    end

                    % remove predictors with zero variance 
                    V = round(var(Ztrial_2D,0,2),6);
                    isUsedPred = V>0;

                    Mdl = fitcdiscr(Ztrial_2D(isUsedPred,:)',isClass,'discrimtype','diaglinear','CrossVal','on'); % 10-fold cross-validation
    %                 Mdl = fitcsvm(Ztrial_2D(isUsedPred,:)',isClass,'CrossVal','on');% 10-fold cross-validation
                    L =  kfoldLoss(Mdl); % Estimate the cross-validated classification error-  the classification errors averaged over all folds

    %                 Random forest
    %                 Mdl = TreeBagger(NumTrees,Ztrial_2D',isClass,'OOBPrediction','On','Method','classification');
    %                 % Mdl.OOBIndices stores the out-of-bag indices as a matrix of logical values
    %                 L = oobError(Mdl);
    %                 correctRate(shuffle,comp,session) = 1 - L(end);
    %                 figure;
    %                 plot(oobErrorBaggedEnsemble)
    %                 xlabel 'Number of grown trees';
    %                 ylabel 'Out-of-bag classification error';

    %                 Optimization of LDA
    %                 Mdl_opt = fitcdiscr(Ztrial_2D(isUsedPred,:)',isClass,'OptimizeHyperparameters','auto',...
    %                         'HyperparameterOptimizationOptions',...
    %                         struct('AcquisitionFunctionName','expected-improvement-plus', 'ShowPlots',false,'Verbose',0,'Kfold',nCV));    
    %                 L = Mdl_opt.HyperparameterOptimizationResults.MinObjective;    
    %                 cvMdl_opt = crossval(Mdl_opt,'KFold',nCV);
    %                 L_opt =  kfoldLoss(cvMdl_opt);

    %                 no CV
    %                 Mdl = fitcdiscr(Ztrial_2D',isClass,'discrimtype','diaglinear'); 
    %                 Mdl = fitcsvm(Ztrial_2D',isClass);
    %                 L = loss(Mdl,Ztrial_2D',isClass);
    %                 
                    correctRate(shuffle,comp,session) = 1 - L;
                    clear Mdl L 
                end
            end
            clear Xtrial numOfTrials 
        end
        Classification_summary{i,1} = recSites{i};
        Classification_summary{i,2} = squeeze(correctRate(1,:,:)); % true 
        Classification_summary{i,3} = squeeze(correctRate(2,:,:)); % shuffled
        % row: component; column: number of session 
        clear correctRate
    end

    % Save classification accuracy 
%     accuracyPath = fullfile(mainDir,'Avg_across_sessions','stimulus_window_0top1','LDA_accuracy_10CV_true&shuffled_mergedL&Rlicks.mat');
    % accuracyPath = fullfile(mainDir,'Avg_across_sessions','SVM_accuracy_10CV_true&shuffled_mergedL&Rlicks.mat');
%     accuracyPath = fullfile(mainDir,'Avg_across_sessions','pre_stimulus_window_-p1to0',...
%         ['RF_accuracy_true&shuffled_mergedL&Rlicks_numTrees',num2str(NumTrees),'.mat']);
    TimeStart = (binNum-1)*bin + timeWindow(1);
    TimeEnd = binNum*bin + timeWindow(1);
    accuracyPath = fullfile(mainDir,'Avg_across_sessions',['LDA_' num2str(TimeStart) 'to' num2str(TimeEnd) '.mat']);
    save(accuracyPath,'Classification_summary')
    clear Classification_summary
end
%% Histogram of true and shuffled classification accuracy 
mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\CorrectTrials\Trial_by_Trial\';
% load([mainDir 'Avg_across_sessions\pre_stimulus_window_-p1to0\RF_accuracy_true&shuffled_mergedL&Rlicks_numTrees500.mat'])
% load([mainDir 'Avg_across_sessions\pre_stimulus_window_-p1to0\SVM_accuracy_10CV_true&shuffled_mergedL&Rlicks.mat'])
% load([mainDir 'Avg_across_sessions\pre_stimulus_window_-p1to0\LDA_accuracy_10CV_true&shuffled_mergedL&Rlicks.mat'])
load([mainDir 'Avg_across_sessions\stimulus_window_0top1\LDA_accuracy_10CV_true&shuffled_buit-inCV_mergedL&Rlicks.mat'])

% block, stimulus, decision 
% Chang et al., Fig S3
figure('Position', [0 0 800 550])
for i = 1: length(recSites)
    for comp = 1:length(component_names)
        for shuffle=1:2
%             subplot(length(recSites),length(component_names),(i-1)*length(component_names)+comp)
            subplot(length(component_names),length(recSites),(comp-1)*length(recSites)+i)
            ca = Classification_summary{i,1+shuffle}(comp,:);
            h = histogram(ca,[0:0.025:1],'FaceColor',shuffle_colors{shuffle},'FaceAlpha',0.6);hold on;
            % label median
            caMedian = median(ca);
            dataMedian(i,comp,shuffle) = caMedian;
            binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end))/2; % Find bin centers
            [~, index] = min(abs(caMedian - binCenters));
            plot([binCenters(index) binCenters(index)],[6.4 8],...
                '-', 'LineWidth', 1, 'Color',shuffle_colors{shuffle});
            if i==length(recSites) && comp == length(component_names)
                text(0.75,(4.5-1*shuffle),shuffle_names{shuffle},'Color',shuffle_colors{shuffle})
            end
        end
        ylim([0 8]);
        set(gca,'box','off','TickDir','out', 'XTick',0:0.25:1, 'YTick',0:2:8)
        xline(0.5,'--');
        if i==1
            ylabel({component_names{comp} 'Number of sessions'})
        end
        if comp==1
            title(recSite_names{i})
        elseif comp==3
            xlabel('Classification accuracy')
        end
%         if comp==1
%             ylabel({recSite_names{i} 'Number of sessions'})
%                 ylabel({recSites{i} 'Number of mice'})
%         end
%         if i==1
%             title(component_names{comp})
%         elseif i==4
%             xlabel('Classification accuracy')
%         end
    end
end
% sgtitle('0-0.1s')
sgtitle('-0.1-0s')
% Save figure
% caFigPath = fullfile(mainDir,'Avg_across_sessions','SVM_accuracy_10CV_true&shuffled_mergedL&Rlicks_histogram');
% caFigPath = fullfile(mainDir,'Avg_across_sessions','stimulus_window_0top1','LDA_accuracy_10CV_true&shuffled_mergedL&Rlicks_histogram');
caFigPath = fullfile(mainDir,'Avg_across_sessions','pre_stimulus_window_-p1to0','RF_accuracy_10CV_true&shuffled_mergedL&Rlicks_histogram');
% caFigPath = fullfile(mainDir,'Avg_across_mice',['RF_accuracy_true&shuffled_mergedL&Rlicks_numTrees',num2str(NumTrees),'_histogram']);
print(caFigPath,'-dpdf','-painters','-loose');


% block
% Chang et al., Fig 3e
figure('Position', [0 0 400 400])
for i = 1: length(recSites)
    subplot(2,2,i)
    for shuffle=1:2
        ca = Classification_summary{i,1+shuffle}(1,:);
        h = histogram(ca,[0:0.025:1],'FaceColor',shuffle_colors{shuffle},'FaceAlpha',0.6);hold on;
        % label median
        caMedian = median(ca);
        dataMedian(i,shuffle) = caMedian;
        binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end))/2; % Find bin centers
        [~, index] = min(abs(caMedian - binCenters));
        plot([binCenters(index) binCenters(index)],[6.4 8],...
            '-', 'LineWidth', 1, 'Color',shuffle_colors{shuffle});
        if i==2
            text(0.7,(6-0.75*shuffle),shuffle_names{shuffle},'Color',shuffle_colors{shuffle})
        end
    end
    ylim([0 8]);
    set(gca,'box','off','TickDir','out', 'XTick',0:0.25:1, 'YTick',0:2:8)
    xline(0.5,'--');
    title(recSite_names{i})
    if i==1 || i==3
        ylabel('Session counts')
    end
    if i==3 || i==4
        xlabel('Classification accuracy')
    end
end
% Save figure
% caFigPath = fullfile(mainDir,'Avg_across_mice','SVM_accuracy_noCV_true&shuffled_mergedL&Rlicks_histogram');
% caFigPath = fullfile(mainDir,'Avg_across_sessions','LDA_accuracy_10CV_true&shuffled_histogram_block');
% caFigPath = fullfile(mainDir,'Avg_across_mice',['RF_accuracy_true&shuffled_mergedL&Rlicks_numTrees',num2str(NumTrees),'_histogram']);
print(caFigPath,'-dpdf','-painters','-loose');

%% Bootstrapped and bootstrapped + shuffled classification accuracy (session)
% Multilevel bootstrapping: sessions -> neurons -> CV partition -> trials 
% To prevent the same trials shown in train and test datasets after bootstrapping, we first perform cross-validation partition and
% then bootstrap trials within train and test datasets seperately. 
% For each bootstrapping, classification accuracy is avergaed across sessions.
% Chang et al., Fig 3f and S3

for binNum= binStart:binEnd   
    for i = 1: length(recSites)
        u = data_all(strcmp(data_all.recSite, recSites{i}), :);
        for boot=1:nBoot  
            %%%%%%%%% Bootstrapping: sessions %%%%%%%%
            SessionInd = randi(size(u,1),size(u,1),1);
            u_resampled = u(SessionInd,:);

            for session=1:size(u_resampled,1)
%                 Xtrial = u_resampled{session,4}{1,1}(:,:,:,:,initialBin,:);
                Xtrial = u_resampled{session,4}{1,1}(:,:,:,:,binNum,:);
                numOfTrials = u_resampled{session,5}{1,1};

                %%%%%%%%% Bootstrapping: neurons %%%%%%%% 
                N = size(Xtrial,1); 
                NeuronInd = randi(N,N,1);
                Xtrial_resampled_n = Xtrial(NeuronInd,:,:,:,:,:);  

                %%%%%%%%% Normalization %%%%%%%%
                Xtrial_2D = [];
                for b = 1:B
                    for s = 1:S
                        for d = 1:D
                            numOfTrials_tt = numOfTrials(1,b,s,d);
                            Xtrial_tt = squeeze(Xtrial_resampled_n(:,b,s,d,:,1:numOfTrials_tt)); 
                            Xtrial_2D = [Xtrial_2D Xtrial_tt];
                        end
                    end
                end
                X_mean = mean(Xtrial_2D,2);
                X_std = std(Xtrial_2D,0,2);
                Ztrial = (Xtrial_resampled_n - X_mean)./(X_std + 0.01); % adding 0.01 to std to prevent the denominator is zero          

                %%%%%%%%% Cross validation partition: train and test datasets %%%%%%%%
                for tt = 1:length(tts)
                    numOfTrials_tt = numOfTrials(1,tts{tt}(1),tts{tt}(2),tts{tt}(3)); 
                    Ztrial_tt = squeeze(Ztrial(:,tts{tt}(1),tts{tt}(2),tts{tt}(3),:,1:numOfTrials_tt)); 
                    % Run this until each partition has >1 trial 
                    for par = 1:50
                        CVParti_ind = randi(nCV,numOfTrials_tt,1);
                        numofTrials_CV = histcounts(CVParti_ind,1:11);
                        TrialInd = sum(numofTrials_CV == 0) == 0;
                        if TrialInd
                            break
                        end
                    end

                    for cv=1:nCV
                        isTest = CVParti_ind==cv;
                        Ztrial_tt_test = Ztrial_tt(:,isTest); 
                        Ztrial_tt_train = Ztrial_tt(:,~isTest);
                        %%%%%%%%% Bootstrapping: trials within train and test datasets %%%%%%%%
                        testTrial_Ind = randi(size(Ztrial_tt_test,2),size(Ztrial_tt_test,2),1);
                        trainTrial_Ind = randi(size(Ztrial_tt_train,2),size(Ztrial_tt_train,2),1);

                        CVParti_data{1,tt,cv} = Ztrial_tt_train(:,trainTrial_Ind);
                        CVParti_data{2,tt,cv} = Ztrial_tt_test(:,testTrial_Ind);
                    end
                end

                %%%%%%%%% Cross Validation %%%%%%%%
                for shuffle=1:length(shuffle_names) % 1: no shuffling; 2: shuffling
                    for comp = 1:length(component_names)
                        for cv=1:nCV  
                            for j=1:2 % j=1: train; j=2: test
                                % concatenation
                                isClass = [];
                                Ztrial_2D = [];
                                for tt = 1:length(tts)
                                    Ztrial_tt = CVParti_data{j,tt,cv};
                                    Ztrial_2D = [Ztrial_2D Ztrial_tt];
                                    if comp ==3 
                                        if tts{tt}(comp) == 2 % merge left and right licks
                                            isClass_tt = ones(size(Ztrial_tt,2),1);
                                        else
                                            isClass_tt = repmat(tts{tt}(comp),[size(Ztrial_tt,2) 1]);
                                        end
                                    else
                                        isClass_tt = repmat(tts{tt}(comp),[size(Ztrial_tt,2) 1]);
                                    end

                                    isClass = [isClass; isClass_tt];
                                    clear isClass_tt Ztrial_tt
                                end
                                if shuffle ==2 % shuffle labels
                                    newOrder = randperm(length(isClass));
                                    isClass = isClass(newOrder);
                                end
                                switch j
                                    case 1 % train
                                        Mdl = fitcdiscr(Ztrial_2D',isClass,'discrimtype','diaglinear'); 
    %                                     Mdl = fitcsvm(Ztrial_2D',isClass);

                                    case 2 % test
                                        prediction = Mdl.predict(Ztrial_2D');
                                        isCorrect = prediction ==isClass;
                                        correctRate(shuffle,comp,session,cv) = sum(isCorrect)/length(isCorrect);
                                        clear Mdl prediction isCorrect 
                                end
                            end
                        end
                    end
                end
                clear Xtrial numOfTrials CVParti_data
            end
            bootstrapped_ca(1:length(shuffle_names),1:length(component_names),boot) = mean(correctRate,[3,4]);
            clear correctRate
        end
        Classification_summary{i,1} = recSites{i};
        Classification_summary{i,2} = squeeze(bootstrapped_ca(1,:,:)); % bootstrapped 
        Classification_summary{i,3} = squeeze(bootstrapped_ca(2,:,:)); % bootstrapped + shuffled 
        % row: component; column: number of bootstrapping 
    end
    % Save classification accuracy 
    % accuracyPath = fullfile(mainDir,'Avg_across_sessions','stimulus_window_0top1',...
    %     ['LDA_accuracy_10CV_bootstrapped&shuffled_', num2str(nBoot),'nBoot_mergedL&Rlicks.mat']);
    TimeStart = (binNum-1)*bin + timeWindow(1);
    TimeEnd = binNum*bin + timeWindow(1);
    accuracyPath = fullfile(mainDir,'Avg_across_sessions',...
        ['LDA_' num2str(nBoot),'nBoot_' num2str(TimeStart) 'to' num2str(TimeEnd) '.mat']);
    save(accuracyPath,'Classification_summary')
    clear Classification_summary
end
%% Plot 95% CI of bootstrapped and shuffled accuracy 
mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\CorrectTrials\Trial_by_Trial\';
% load([mainDir 'Avg_across_sessions\pre_stimulus_window_-p1to0\LDA_accuracy_10CV_bootstrapped&shuffled_1000nBoot_mergedL&Rlicks'])
load([mainDir 'Avg_across_sessions\stimulus_window_0top1\LDA_accuracy_10CV_bootstrapped&shuffled_1000nBoot_mergedL&Rlicks'])

figure('Position', [0 0 300 1000])
for i = 1: length(recSites)
    for comp = 1:length(component_names)
        subplot(length(component_names),1,comp)
        for shuffle=1:2  
%             accuracy = sort(squeeze(accuracy_all{i,r+3}(comp,:,:)),2); 
            % bootstrapped accuracy: 4th column; bootstrapped shuffled accuracy: 5th column
            accuracy = sort(squeeze(Classification_summary{i,1+shuffle}(comp,:)));
            accuracy_mean = mean(accuracy);
            accuracy_up95CI = accuracy(round(0.975*nBoot));
            accuracy_low95CI = accuracy(round(0.025*nBoot));
            dataAccuracy(i,comp,shuffle,1) =  accuracy_up95CI;
             dataAccuracy(i,comp,shuffle,2) =  accuracy_low95CI;
            errorbar(i+shuffle*0.2-0.3,accuracy_mean,(accuracy_mean-accuracy_low95CI),(accuracy_up95CI-accuracy_mean),'o',...
                'Color', shuffle_colors{shuffle},'MarkerEdgeColor', shuffle_colors{shuffle},...
                'MarkerFaceColor',shuffle_colors{shuffle},'CapSize', 6, 'LineWidth', 0.5); hold on;
            if i==length(recSites) && comp == 1
                text(3.5,(0.4-0.05*shuffle),shuffle_names{shuffle},'Color',shuffle_colors{shuffle})
            end
        end
        set(gca,'box','off','TickDir','out')
        xlim([0.5 4.5]);ylim([0.2 1]);
        xticks([1:4]);xticklabels(recSite_names);
        ylabel('Classification accuracy');
        yline(0.5,'--');
        title(component_names{comp})
    end
end
% Save figure
caFigPath = fullfile(mainDir,'Avg_across_sessions','pre_stimulus_window_-p1to0',...
    ['LDA_accuracy_bootstrapped&shuffled_', num2str(nBoot),'nBoot_mergedL&Rlicks']);
% caFigPath = fullfile(mainDir,'Avg_across_sessions','CVpartition_after_bootstrapping',['RF_accuracy_10CV_bootstrapped&shuffled2_built-inCV', num2str(nBoot),'nBoot_mergedL&Rlicks']);
print(caFigPath,'-dpdf','-painters','-loose');

%% Example session
% Chang et al., Fig 3d
% The fitcdiscr (matlab nuilt-in function for LDA) function uses Bayesian approach not dimension reduction, 
% so the projeciton matrix is not in the outputs of fitcdiscr. This approach focuses on classificaiton not dimension reduction.
% To visualize the projections, dimension reduction is applied. 

% Reference
% https://stats.stackexchange.com/questions/111421/reproduce-linear-discriminant-analysis-projection-plot
% https://www.mathworks.com/help/stats/classificationdiscriminant-class.html
% https://sebastianraschka.com/Articles/2014_python_lda.html

% Explained variance
% https://stats.stackexchange.com/questions/67342/proportion-of-explained-variance-in-pca-and-lda
% https://stats.stackexchange.com/questions/48786/algebra-of-lda-fisher-discrimination-power-of-a-variable-and-linear-discriminan/48859#48859
% Here, the discriminability (signal-to-noise ratio) of latent variables is more important than their captured variances 
% since we focused on their roles in discriminating classes not in explaining variance of data.

% colors = {[0.4 0.85 0.1] [1 0.4 0.16]};
colors = {[0.05 0.5 0.1] [0.9 0.8 0.2]};

for i = 1: length(recSites)
    recSite = recSites{i};
    u = table2cell(data_all(strcmp(data_all.recSite, recSite), :));
    for session=1:size(u,1)
        Xtrial = u{session,4}(:,:,:,:,initialBin,:);
        numOfTrials = u{session,5};
        N = size(Xtrial,1); % number of neurons
        %%%%%%%%% Normalization %%%%%%%%
        Xtrial_2D = [];
        for b = 1:B
            for s = 1:S
                for d = 1:D
                    numOfTrials_tt = numOfTrials(1,b,s,d);
                    Xtrial_tt = squeeze(Xtrial(:,b,s,d,:,1:numOfTrials_tt)); 
                    Xtrial_2D = [Xtrial_2D Xtrial_tt];
                end
            end
        end
        X_mean = mean(Xtrial_2D,2);
        X_std = std(Xtrial_2D,0,2);
        Ztrial = (Xtrial - X_mean)./(X_std + 0.01); % adding 0.01 to std to prevent the denominator is zero  
        %%%%%%%%% Concatenation %%%%%%%%
        isClass = [];
        Ztrial_2D = [];
        for tt = 1:length(tts)
            numOfTrials_tt = numOfTrials(1,tts{tt}(1),tts{tt}(2),tts{tt}(3)); 
            Ztrial_tt = squeeze(Ztrial(:,tts{tt}(1),tts{tt}(2),tts{tt}(3),1:numOfTrials_tt));
            Ztrial_2D = [Ztrial_2D Ztrial_tt];
            isClass_tt = repmat(tts{tt}(1),[numOfTrials_tt 1]);
            isClass = [isClass; isClass_tt];
            clear isClass_tt Ztrial_tt
        end 
        % remove predictors with zero variance 
        J = round(var(Ztrial_2D,0,2),6);
        isUsedPred = J>0;
        Ztrial_2D = Ztrial_2D(isUsedPred,:)';
        
        %%%%%%%%% Built-in LDA function %%%%%%%%
        Mdl_CV = fitcdiscr(Ztrial_2D,isClass,'discrimtype','diaglinear','CrossVal','on'); % 10-fold cross-validation
        L =  kfoldLoss(Mdl_CV);
        correctRate = 1-L;
        Mdl = fitcdiscr(Ztrial_2D,isClass,'discrimtype','linear'); 
        Sw = Mdl.Sigma;
        Sb = Mdl.BetweenSigma;
        
%         %%%%%%%%% Step-by-step LDA calculation %%%%%%%%
%         for g=1:max(isClass) % loop over groups
%             mus(g,:) = mean(Ztrial_2D(isClass==g,:)); % # class means
%             Ng(g) = length(find(isClass==g)); % # number of points per group
%         end
%         mu = mean(Ztrial_2D); % overall mean
%         Sw = zeros(size(Ztrial_2D,2)); % # within-class scatter matrix
%         Sb = zeros(size(Ztrial_2D,2)); % # between-class scatter matrix
%         for g=1:max(isClass)
%             Xg = bsxfun(@minus, Ztrial_2D(isClass==g,:), mus(g,:)); % # centred group data
%             Sw = Sw + transpose(Xg)*Xg;
%             Sb = Sb + Ng(g)*(transpose(mus(g,:) - mu)*(mus(g,:) - mu));
%         end
%         St = transpose(bsxfun(@minus,Ztrial_2D,mu)) * bsxfun(@minus,Ztrial_2D,mu); % # total scatter matrix
%         assert(sum(sum((St-Sw-Sb).^2)) < 1e-10, 'Error: Sw + Sb ~= St')
        
        % LDA 
        [V,E] = eig(Sw\Sb); % Note: if the input matrix is not symmetric, it will return complex number
        [e, ind] = sort(diag(E),'descend'); 
        Es = E(ind,ind); Vs = V(:,ind);
        Z =  Ztrial_2D*Vs;
        
        % session info
        mouseName = u{session,1};
        seshDate = u{session,2};
        %%%%%%%%% Plot %%%%%%%%
        if isreal(Z(:,1)) && isreal(Z(:,2))
            % Projections onto the LD1 and LD2 coordinates
            figure('Position',[0 0 800 800])
            subplot('Position',[0.1 0.1 0.6 0.6])
            for g=1:max(isClass)
                scatter(Z(isClass==g,1),Z(isClass==g,2), 'MarkerEdgeColor',colors{g}); hold on;
            end
            set(gca, 'box','off','TickDir','out')
            xlabel('LD1','FontSize',10);
            ylabel('LD2','FontSize',10);
            xl = xlim; yl = ylim;
            xbin = round(diff(xl)/20,1);
            ybin = round(diff(yl)/20,1);
            % Histogram onto the LD1 axis 
            subplot('Position',[0.1 0.75 0.6 0.15])
            for g=1:max(isClass)
                histogram(Z(isClass==g,1), 'Normalization', 'probability', 'binWidth', xbin, 'FaceColor', colors{g}); hold on;
            end
            set(gca, 'box','off','TickDir','out')
            xlim(xl);
            xticklabels([]); 
            ylabel('Proportion of trials','FontSize',10);

            % Histogram onto the LD1 axis 
            subplot('Position',[0.75 0.1 0.15 0.6])
            for g=1:max(isClass)
                histogram(Z(isClass==g,2), 'Normalization', 'probability', 'binWidth', ybin, 'FaceColor', colors{g}, ...
                    'Orientation','horizontal'); hold on;
            end
            set(gca, 'box','off','TickDir','out')
            ylim(yl);
            yticklabels([]); 
            xlabel('Proportion of trials','FontSize',10);
            
        elseif isreal(Z(:,1))
            % Histogram onto the LD1 axis
            figure('Position',[0 0 400 400])
            for g=1:max(isClass)
                histogram(Z(isClass==g,1), 'Normalization', 'probability', 'binWidth', 0.4, 'FaceColor', colors{g}); hold on;
            end
            set(gca, 'box','off','TickDir','out')
            xlabel('LD1','FontSize',10);
            ylabel('Proportion of trials','FontSize',10);
        end
        
        sgtitle({[recSite '\_' mouseName '\_' seshDate]...
            ['accuracy = ' num2str(round(correctRate,2)) '; n = ' num2str(N)]})
        
        % Save figure
        ldaFigPath = fullfile(mainDir,'Example_sessions', recSite,[mouseName '_' seshDate]);
        print(ldaFigPath,'-dpdf','-painters','-loose');    
        clear Xtrial numOfTrials Sw Sb St mus Ng
    end  
    close all
end

%% Plot histogram  across time
mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\CorrectTrials\Trial_by_trial';
bin = 0.1;
timeWindow = [-1 2.5];
analysisWindow = [-1 0.3];
binStart =round((analysisWindow(1)-timeWindow(1))/bin+1);
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
shuffle_names = {'true' 'shuffled'};
shuffle_colors = {[0 1 1] [0.5 0.5 0.5]};

for binNum= binStart:binEnd
    TimeStart = (binNum-1)*bin + timeWindow(1);
    TimeEnd = binNum*bin + timeWindow(1);
    accuracyPath = fullfile(mainDir,'Avg_across_sessions', 'LDA_-1top3', ['LDA_' num2str(TimeStart) 'to' num2str(TimeEnd) '.mat']);
    load(accuracyPath)
    figure('Position', [0 0 800 200])
    for i = 1: length(recSites)
        subplot(1,4,i)
        for shuffle=1:2
            ca = Classification_summary{i,1+shuffle}(1,:);
            h = histogram(ca,[0:0.025:1],'FaceColor',shuffle_colors{shuffle},'FaceAlpha',0.6);hold on;
            % label median
            caMedian = median(ca);
            dataMedian(i,shuffle) = caMedian;
            binCenters = (h.BinEdges(1:end-1) + h.BinEdges(2:end))/2; % Find bin centers
            [~, index] = min(abs(caMedian - binCenters));
            plot([binCenters(index) binCenters(index)],[6.4 8],...
                '-', 'LineWidth', 1, 'Color',shuffle_colors{shuffle});
            if i==4
                text(0.7,(7-0.75*shuffle),shuffle_names{shuffle},'Color',shuffle_colors{shuffle})
            end
        end
        ylim([0 8]);
        set(gca,'box','off','TickDir','out', 'XTick',0:0.25:1, 'YTick',0:2:8)
        xline(0.5,'--');
        title(recSite_names{i})
        if i==1 
            ylabel('Session counts')
        end
        xlabel('Classification accuracy')
    end
    sgtitle([num2str(TimeStart) ' to ' num2str(TimeEnd)])
    % Save figure
    caFigPath = fullfile(mainDir,'Avg_across_sessions','LDA_-1top3',['LDA_histogram_' num2str(TimeStart*1000) 'to' num2str(TimeEnd*1000)]);
    print(caFigPath,'-dpdf','-painters','-loose');
end

%% Plot bootstrapped data across time

mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\CorrectTrials\Trial_by_Trial\';
bin = 0.1;
timeWindow = [-1 2.5];
analysisWindow = [-1 0.3];
binStart =round((analysisWindow(1)-timeWindow(1))/bin+1);
binEnd = round((analysisWindow(2)-timeWindow(1))/bin);
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1', 'S2', 'MM', 'ALM'};
shuffle_names = {'true' 'shuffled'};
shuffle_colors = {[0 1 1] [0.5 0.5 0.5]};
component_names = {'Block' 'Stimulus' 'Decision'};
nBoot = 100;
for binNum= binStart:binEnd
    TimeStart = (binNum-1)*bin + timeWindow(1);
    TimeEnd = binNum*bin + timeWindow(1);
    accuracyPath = fullfile(mainDir,'Avg_across_sessions', 'LDA_100nBoot_-1top3', ['LDA_100nBoot_' num2str(TimeStart) 'to' num2str(TimeEnd) '.mat']);
    load(accuracyPath)
    figure('Position', [0 0 300 1000])
    for i = 1: length(recSites)
        for comp = 1:length(component_names)
            subplot(length(component_names),1,comp)
            for shuffle=1:2  
    %             accuracy = sort(squeeze(accuracy_all{i,r+3}(comp,:,:)),2); 
                % bootstrapped accuracy: 4th column; bootstrapped shuffled accuracy: 5th column
                accuracy = sort(squeeze(Classification_summary{i,1+shuffle}(comp,:)));
                accuracy_mean = mean(accuracy);
                accuracy_up95CI = accuracy(round(0.975*nBoot));
                accuracy_low95CI = accuracy(round(0.025*nBoot));
                dataAccuracy(i,comp,shuffle,1) =  accuracy_up95CI;
                 dataAccuracy(i,comp,shuffle,2) =  accuracy_low95CI;
                errorbar(i+shuffle*0.2-0.3,accuracy_mean,(accuracy_mean-accuracy_low95CI),(accuracy_up95CI-accuracy_mean),'o',...
                    'Color', shuffle_colors{shuffle},'MarkerEdgeColor', shuffle_colors{shuffle},...
                    'MarkerFaceColor',shuffle_colors{shuffle},'CapSize', 6, 'LineWidth', 0.5); hold on;
                if i==length(recSites) && comp == 1
                    text(3.5,(0.4-0.05*shuffle),shuffle_names{shuffle},'Color',shuffle_colors{shuffle})
                end
                
                ca_mean(i,comp,shuffle,binNum) = accuracy_mean;
                ca_up95CI(i,comp,shuffle,binNum) = accuracy_up95CI; 
                ca_low95CI(i,comp,shuffle,binNum) = accuracy_low95CI;
            end
            set(gca,'box','off','TickDir','out')
            xlim([0.5 4.5]);ylim([0.2 1]);
            xticks([1:4]);xticklabels(recSite_names);
            ylabel('Classification accuracy');
            yline(0.5,'--');
            title(component_names{comp})
        end
    end
    sgtitle([num2str(TimeStart) ' to ' num2str(TimeEnd)])
    % Save figure
    caFigPath = fullfile(mainDir,'Avg_across_sessions','LDA_100nBoot_-1top3',...
        ['LDA_100nBoot_' num2str(TimeStart*1000) 'to' num2str(TimeEnd*1000)]);
    print(caFigPath,'-dpdf','-painters','-loose');
end

time = analysisWindow(1)+bin:bin:analysisWindow(2);
time_bw = [time, fliplr(time)];
stimWindow = [0 0.15];
plotWindow = analysisWindow;
for comp = 1:length(component_names)
    figure('Position', [0 0 600 600])
    for i = 1: length(recSites)
        subplot(2,2,i)
        for shuffle=1:2
            plot(time, squeeze(ca_mean(i,comp,shuffle,:)),'Color',shuffle_colors{shuffle}); hold on;
            fill(time_bw,[squeeze(ca_low95CI(i,comp,shuffle,:))', fliplr(squeeze(ca_up95CI(i,comp,shuffle,:))')],...
                shuffle_colors{shuffle},'FaceAlpha', 0.3, 'Linestyle', 'none');
            plot(stimWindow, [1 1],'k','Linewidth',2)
            set(gca, 'box','off','TickDir','out')
            xlim(plotWindow);ylim([0.4 1]);
            ylabel('Classification accuracy');
            xlabel('Time from stimulus onset (s)')
            title(recSite_names{i})
        end
    end
    sgtitle(component_names{comp})
    % Save figure
    caFigPath = fullfile(mainDir,'Avg_across_sessions','LDA_100nBoot_-1top3',...
        ['LDA_100nBoot_' component_names{comp}]);
    print(caFigPath,'-dpdf','-painters','-loose');
end


        