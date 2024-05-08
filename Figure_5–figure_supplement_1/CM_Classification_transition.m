% Classification of initial states of transition trials
% Trial-by-trial analysis
% Initial state: -100 ~ 0 ms activity (100 ms before stimulus onsets)
% We define transition as the time period from the first trial to the first hit trial. 

% Are the initial states of transition trials closer to the initial states of previous block or current block?


% True prediction
% 1. Create the training (90%) and testing (10%) sets for correct trials after transition
% 2. Train the classifier
% 3. Predict the class of transition trials and the testing set

% load data array
load('E:\Data_array\AllTrials_100bin_noSmooth\data_array.mat')
% Import setting 
mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\Transition\';
B=2;
S=2;
D=3;
bin = 0.1;
initialBin = (0+1)*(1/bin);
recSites = {'left S1', 'left S2', 'left wM2','left ALM'};

% % tBlock vs. vBlock
train_tts = {[1 1 1] [1 2 3] [2 1 3] [2 2 2]}; % trial types for training: TTT (tHit), TVN (vCR), VTN (tCR), VVV (vHit) 
comp = 1; % 1: block, 2: stimulus, 3: decision
Class1_name = 'respond-to-touch'; % name of class 1 
Class2_name = 'respond-to-light'; % name of class 2
% shuffle_names = {'true' 'shuffled'};
nBoot = 1000;
nCV = 5;
NumTrees = 500;
%% True prediction (session)
for i = 1: length(recSites)
    u = table2cell(data_all(strcmp(data_all.recSite, recSites{i}), :));
    for session=1:size(u,1)
        Xtrial = u{session,4}(:,:,:,:,initialBin,:);
        numOfTrials = u{session,5};
        early_transition = u{session,6};
        late_transition = u{session,7};
        not_transition = u{session,8};
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
        
        %%%%%%%%% Training/testing sets %%%%%%%% 
        % Concatenation 
        isClass_train_2D = [];
        Ztrial_train_2D = [];
        for train = 1:length(train_tts)
            numOfTrials_tt = numOfTrials(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3)); 
%             early_transition_tt = squeeze(early_transition(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),1:numOfTrials_tt));
%             late_transition_tt = squeeze(late_transition(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),1:numOfTrials_tt));
            not_transition_tt = logical(squeeze(not_transition(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),1:numOfTrials_tt)));
            Ztrial_tt = squeeze(Ztrial(:,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),:,not_transition_tt));
%             transition_tt = early_transition_tt | late_transition_tt;
%             Ztrial_tt = squeeze(Ztrial(:,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),:,~transition_tt));
            
            newInd = randperm(size(Ztrial_tt,2));
            numofTesting = round(size(Ztrial_tt,2)/nCV);
            
            Ind_test = newInd <= numofTesting;
            Ind_train = newInd > numofTesting;
            Ztrial_tt_test = Ztrial_tt(:,Ind_test);
            Ztrial_tt_train = Ztrial_tt(:,Ind_train);
            Ztrial_train_2D = [Ztrial_train_2D Ztrial_tt_train];
            Ztrial_train{train} = Ztrial_tt_train;
            Ztrial_test{train} = Ztrial_tt_test;
            
            isClass_tt_train = repmat(train_tts{train}(comp),[sum(Ind_train) 1]);
            isClass_train_2D = [isClass_train_2D; isClass_tt_train];
            isClass_train{train} = isClass_tt_train;
            isClass_tt_test = repmat(train_tts{train}(comp),[sum(Ind_test) 1]);
            isClass_test{train} = isClass_tt_test;
            
            clear Ztrial_tt Ztrial_tt_test Ztrial_tt_train isClass_tt_train isClass_tt_test numOfTrials_tt not_transition_tt...
%                 early_transition_tt late_transition_tt transition_tt
        end
        % remove predictors with zero variance 
        V = round(var(Ztrial_train_2D,0,2),6);
        isUsedPred = V>0;
        %%%%%%%%% Training %%%%%%%%
%         Mdl = fitcdiscr(Ztrial_train_2D(isUsedPred,:)',isClass_train_2D,'discrimtype','linear');
%         mdl = fitcdiscr(Ztrial_2D_train',isClass_train,'discrimtype','diagQuadratic');
%         Mdl = fitcsvm(Ztrial_train_2D(isUsedPred,:)',isClass_train_2D);
%         CompactMdl = compact(Mdl);
%         % Estimate the optimal score-to-posterior-probability transformation function
%         CompactMdl = fitPosterior(CompactMdl,Ztrial_2D_train(isUsedPred,:)',isClass_train);
%           Random forest
%         Mdl = TreeBagger(NumTrees,Ztrial_2D_train(isUsedPred,:)',isClass_train,'Method','classification');
%       Optimization of LDA
        Mdl = fitcdiscr(Ztrial_train_2D(isUsedPred,:)',isClass_train_2D,'OptimizeHyperparameters','auto',...
                'HyperparameterOptimizationOptions',...
                struct('AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0,'Kfold',nCV)); 

%         % Optimization of SVM
%         Mdl = fitcsvm(Ztrial_train_2D(isUsedPred,:)',isClass_train_2D,'OptimizeHyperparameters','auto',...
%             'HyperparameterOptimizationOptions',...
%             struct('AcquisitionFunctionName','expected-improvement-plus','Kfold',nCV));
%         CompactMdl = compact(Mdl);
%         CompactMdl = fitPosterior(CompactMdl,Ztrial_train_2D(isUsedPred,:)',isClass_train_2D);
       
        %%%%%%%%% Testing for correct trails %%%%%%%%
        for train = 1:length(train_tts)
            Ztrial_tt_test = Ztrial_test{train};
            [labels, postProbs] = Mdl.predict(Ztrial_tt_test(isUsedPred,:)');
%             [labels, postProbs] = CompactMdl.predict(Ztrial_tt_test(isUsedPred,:)');
            isClass1 = labels ==1;
%             % Random forest
%             isClass1 = cell2mat(labels) == '1';
            test_label(train,session) =  sum(isClass1)/length(isClass1); % tBlock label
            test_probability(train,session) =  mean(postProbs(:,1)); % tBlock probability
            clear Ztrial_tt_test 
        end
        
         %%%%%%%%% Prediction for transition trials %%%%%%%%
        % labels: predicted class labels
        % postProbs: predicted posterior probabilities- the probability that an observation belongs in a particular class
        for g=1:2
            switch g 
                case 1
                    transition = early_transition;
                case 2
                    transition = late_transition;
            end

            for b=1:B 
                Ztrial_block = squeeze(Ztrial(:,b,:,:,:,:));
                transition_block = squeeze(transition(:,b,:,:,:));
                Ztrial_pred = reshape(Ztrial_block(logical(transition_block)),N,[]);
                [labels, postProbs] = Mdl.predict(Ztrial_pred(isUsedPred,:)');
%                 [labels, postProbs] = CompactMdl.predict(Ztrial_pred(isUsedPred,:)');
                isClass1 = labels ==1;
    %             % Random forest
    %             isClass1 = cell2mat(labels) == '1';
                prediction_label((g-1)*B+b,session) =  sum(isClass1)/length(isClass1); 
                prediction_probability((g-1)*B+b,session) =  mean(postProbs(:,1)); 
                clear Ztrial_block transition_block Ztrial_pred
            end
        end
        
        %%%%%%%%% Positive control for correct trails %%%%%%%%
        for train = 1:length(train_tts)
            Ztrial_tt_train = Ztrial_train{train};
            [labels, postProbs] = Mdl.predict(Ztrial_tt_train(isUsedPred,:)');
%             [labels, postProbs] = CompactMdl.predict(Ztrial_tt_train(isUsedPred,:)');
            isClass1 = labels ==1;
%             % Random forest
%             isClass1 = cell2mat(labels) == '1';
            train_label(train,session) =  sum(isClass1)/length(isClass1); % tBlock label
            train_probability(train,session) =  mean(postProbs(:,1)); % tBlock probability
            clear Ztrial_tt_train 
        end
        clear Xtrial numOfTrials early_transition late_transition transition
    end
    Prediction_summary{i,1} = recSites{i};
    Prediction_summary{i,2} = test_label; % row: TTT (tHit), TVN (vCR), VTN (tCR), VVV (vHit) ; column: sessions
    Prediction_summary{i,3} = test_probability;
    Prediction_summary{i,4} = prediction_label;     % row: early trials in tBlock, early trials in vBlock, late trials in tBlock, late trials in vBlock
    Prediction_summary{i,5} = prediction_probability; 
    Prediction_summary{i,6} = train_label; 
    Prediction_summary{i,7} = train_probability; 
    clear test_label test_probability prediction_label prediction_probability train_label train_probability
end
predPath = fullfile(mainDir,['LDA_true_prediction_' num2str(nCV) 'CV_optimization.mat']);
% predPath = fullfile(mainDir,'Avg_across_sessions',['SVM_true_prediction_' num2str(nCV) 'CV_allTrials_optimization_-100-0ms.mat']);
% predPath = fullfile(mainDir,'Avg_across_sessions',['RF_true_prediction_numTrees',num2str(NumTrees),'_allTrials_-100-0ms.mat']);
save(predPath,'Prediction_summary')        

%% Bootstrapped prediction (session)
% Multilevel bootstrapping: sessions -> neurons -> trials 
% For each bootstrapping, classification accuracy is avergaed across sessions.
        
for i = 1: length(recSites)
    rng(5)
    u = data_all(strcmp(data_all.recSite, recSites{i}), :);
    for boot=1:nBoot  
        %%%%%%%%% Bootstrapping: sessions %%%%%%%%
        SessionInd = randi(size(u,1),size(u,1),1);
        u_resampled = u(SessionInd,:);
        numOfSessions = size(u_resampled,1);
        for session=1:numOfSessions
            Xtrial = u_resampled{session,4}{1,1}(:,:,:,:,initialBin,:);
            numOfTrials = u_resampled{session,5}{1,1};
            early_transition = u_resampled{session,6}{1,1};
            late_transition = u_resampled{session,7}{1,1};
            not_transition = u_resampled{session,8}{1,1};
            %%%%%%%%% Bootstrapping: neurons %%%%%%%% 
            N = size(Xtrial,1); 
            NeuronInd = randi(N,N,1);
            Xtrial_resampled_n = Xtrial(NeuronInd,:,:,:,:,:);  

            %%%%%%%%% Normalization %%%%%%%%
            Xtrial_2D = [];
            for b =1:B
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
                            
            %%%%%%%%% Training/testing sets %%%%%%%% 
            % Concatenation 
            isClass_train_2D = [];
            Ztrial_train_2D = [];
            for train = 1:length(train_tts)
                numOfTrials_tt = numOfTrials(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3));
                not_transition_tt = logical(squeeze(not_transition(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),1:numOfTrials_tt)));
                Ztrial_tt = squeeze(Ztrial(:,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),:,not_transition_tt));                
%                 early_transition_tt = squeeze(early_transition(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),1:numOfTrials_tt));
%                 late_transition_tt = squeeze(late_transition(1,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),1:numOfTrials_tt));
%                 transition_tt = early_transition_tt | late_transition_tt;
%                 Ztrial_tt = squeeze(Ztrial(:,train_tts{train}(1),train_tts{train}(2),train_tts{train}(3),:,~transition_tt));
                
                newInd = randperm(size(Ztrial_tt,2));
                numofTesting = round(size(Ztrial_tt,2)/nCV);
                
                Ind_test = newInd <= numofTesting;
                Ind_train = newInd > numofTesting;
                Ztrial_tt_test = Ztrial_tt(:,Ind_test);
                Ztrial_tt_train = Ztrial_tt(:,Ind_train);
                
                %%%%%%%%% Bootstrapping: trials %%%%%%%%
                Ind_test_resampled = randi(sum(Ind_test),sum(Ind_test),1);
                Ztrial_test_resampled = Ztrial_tt_test(:,Ind_test_resampled);
                Ind_train_resampled = randi(sum(Ind_train),sum(Ind_train),1);
                Ztrial_train_resampled = Ztrial_tt_train(:,Ind_train_resampled);

                Ztrial_train_2D = [Ztrial_train_2D Ztrial_train_resampled];
                Ztrial_train{train} = Ztrial_train_resampled;
                Ztrial_test{train} = Ztrial_test_resampled;
                
                isClass_tt_train = repmat(train_tts{train}(comp),[sum(Ind_train) 1]);
                isClass_train_2D = [isClass_train_2D; isClass_tt_train];
                isClass_train{train} = isClass_tt_train;
                isClass_tt_test = repmat(train_tts{train}(comp),[sum(Ind_test) 1]);
                isClass_test{train} = isClass_tt_test;

                clear Ztrial_tt Ztrial_tt_test Ztrial_tt_train isClass_tt_train isClass_tt_test numOfTrials_tt...
                    early_transition_tt late_transition_tt transition_tt
            end   
            % remove predictors with zero variance 
            V = round(var(Ztrial_train_2D,0,2),6);
            isUsedPred = V>0;            
            %%%%%%%%% Training %%%%%%%%
            Mdl = fitcdiscr(Ztrial_train_2D(isUsedPred,:)',isClass_train_2D,'discrimtype','linear');
%         mdl = fitcdiscr(Ztrial_2D_train',isClass_train,'discrimtype','diagQuadratic');
%             Mdl = fitcsvm(Ztrial_2D_train(isUsedPred,:)',isClass_train);
%             CompactMdl = compact(Mdl);
%             % Estimate the optimal score-to-posterior-probability transformation function
%             CompactMdl = fitPosterior(CompactMdl,Ztrial_2D_train(isUsedPred,:)',isClass_train);

%             Optimization of LDA
%             Mdl = fitcdiscr(Ztrial_train_2D(isUsedPred,:)',isClass_train_2D,'OptimizeHyperparameters','auto',...
%                     'HyperparameterOptimizationOptions',...
%                     struct('AcquisitionFunctionName','expected-improvement-plus','ShowPlots',false,'Verbose',0,'Kfold',nCV));
                
            %%%%%%%%% Testing for correct trails %%%%%%%%
            for train = 1:length(train_tts)
                Ztrial_tt_test = Ztrial_test{train};
                [labels, postProbs] = Mdl.predict(Ztrial_tt_test(isUsedPred,:)');
%                 [labels, postProbs] = CompactMdl.predict(Ztrial_tt_test(isUsedPred,:)');
                isClass1 = labels ==1;
    %             % Random forest
    %             isClass1 = cell2mat(labels) == '1';
                test_label(train,session) =  sum(isClass1)/length(isClass1); % tBlock label
                test_probability(train,session) =  mean(postProbs(:,1)); % tBlock probability
                clear Ztrial_tt_test 
            end

            %%%%%%%%% Prediction for transition trials %%%%%%%%
            % labels: predicted class labels
            % postProbs: predicted posterior probabilities- the probability that an observation belongs in a particular class
            for g=1:2
                switch g 
                    case 1
                        transition = early_transition;
                    case 2
                        transition = late_transition;
                end
                for b = 1:B
                    Ztrial_block = squeeze(Ztrial(:,b,:,:,:,:));
                    transition_block = squeeze(transition(:,b,:,:,:));
                    Ztrial_pred = reshape(Ztrial_block(logical(transition_block)),N,[]);
                    if ~isempty(Ztrial_pred)
                        %%%%%%%%% Bootstrapping: trials %%%%%%%%
                        Ind_pred_resampled = randi(size(Ztrial_pred,2),size(Ztrial_pred,2),1);
                        Ztrial_pred_resampled = Ztrial_pred(:,Ind_pred_resampled);

                        [labels, postProbs] = Mdl.predict(Ztrial_pred_resampled(isUsedPred,:)');
        %                 [labels, postProbs] = CompactMdl.predict(Ztrial_pred(isUsedPred,:)');
                        isClass1 = labels ==1;
                        prediction_label((g-1)*B+b,session) =  sum(isClass1)/length(isClass1); 
                        prediction_probability((g-1)*B+b,session) =  mean(postProbs(:,1)); 
                    else
                        prediction_label((g-1)*B+b,session) =  NaN;
                        prediction_probability((g-1)*B+b,session) = NaN;
                    end
                    clear Ztrial_block transition_block Ztrial_pred
                end
            end
            %%%%%%%%% Positive control for correct trails %%%%%%%%
            for train = 1:length(train_tts)
                Ztrial_tt_train = Ztrial_train{train};
                [labels, postProbs] = Mdl.predict(Ztrial_tt_train(isUsedPred,:)');
%                 [labels, postProbs] = CompactMdl.predict(Ztrial_tt_test(isUsedPred,:)');
                isClass1 = labels ==1;
                train_label(train,session) =  sum(isClass1)/length(isClass1); % tBlock label
                train_probability(train,session) =  mean(postProbs(:,1)); % tBlock probability
                clear Ztrial_tt_train 
            end
            
            clear Xtrial numOfTrials 
        end
        
        bootstrapped_test_label(1:length(train_tts),1:numOfSessions,boot) = test_label;
        bootstrapped_test_prob(1:length(train_tts),1:numOfSessions,boot) = test_probability;
        bootstrapped_prediction_label(1:B*2, 1:numOfSessions,boot) = prediction_label;
        bootstrapped_prediction_prob(1:B*2,1:numOfSessions,boot) = prediction_probability;
        bootstrapped_train_label(1:length(train_tts),1:numOfSessions,boot) = train_label;
        bootstrapped_train_prob(1:length(train_tts),1:numOfSessions,boot) = train_probability;

        clear test_label test_probability prediction_label prediction_probability train_label train_probability 
    end
    
    Prediction_summary{i,1} = recSites{i};
    Prediction_summary{i,2} = bootstrapped_test_label;
    Prediction_summary{i,3} = bootstrapped_test_prob;
    Prediction_summary{i,4} = bootstrapped_prediction_label; 
    Prediction_summary{i,5} = bootstrapped_prediction_prob; 
    Prediction_summary{i,6} = bootstrapped_train_label; 
    Prediction_summary{i,7} = bootstrapped_train_prob; 
    clear bootstrapped_test_label bootstrapped_test_prob bootstrapped_prediction_label bootstrapped_prediction_prob...
        bootstrapped_train_label trapped_train_prob
end
% Save classification accuracy 
predPath = fullfile(mainDir,['LDA_bootstrapped_prediction_', num2str(nCV),'CV_', num2str(nBoot),'nBoot.mat']);
save(predPath,'Prediction_summary')

%% Kendall's Tau Coefficient and Plotting
% rank correlation
% one tail
% Kendall's tau
% tBlock probabilities of each session are used to calculate Kendall's tau. (size: number of sessions)
% Chang et al., Fig 4c-d
% load data
load('E:\CM_NeuralActivity_Analysis\Classification\Transition\LDA_true_prediction_10CV_optimization') % true prediction
trueData = Prediction_summary;
load('E:\CM_NeuralActivity_Analysis\Classification\Transition\LDA_bootstrapped_prediction_10CV_1000nBoot') % bootstrapped prediction
bootstrpData = Prediction_summary;
mainDir = 'E:\CM_NeuralActivity_Analysis\Classification\Transition';

recSites = {'left S1', 'left S2', 'left wM2','left ALM'};
recSite_names = {'S1' 'S2' 'MM' 'ALM'};
B=2;
nCV=10;
nBoot = 1000;
Class1_name = 'Respond-to-touch'; % name of class 1 
Class2_name = 'Respond-to-light'; % name of class 2
Colors_block = {[0.05 0.5 0.1] [0.9 0.8 0.2]}; % respond-to-touch and respond-to-light block
% Colors_block = {[0.4 0.85 0.1] [1 0.4 0.15]}; % respond-to-touch and respond-to-light block
Colors_transition = {[0.7 0.65 0.75] [0.5 0.2 0.55]}; % early and late
pred_names = {'Light to touch' 'Touch to light'};
xtickNames = {'Pre-' 'Early' 'Late' 'Post'};
for i=1:length(recSites)
    figure('Position', [0 0 600 600])
    %%% True prediction
%     idx = sum(isnan(trueData{i,5}),1)==0; % Remove sessions without early or late trials
%     prediction_probability = trueData{i,5}(:,idx); % row: early trials in tBlock, early trials in vBlock, late trials in tBlock, late trials in vBlock
%     test_tBlock = mean(trueData{i,3}(1:2,idx),1);
%     test_vBlock = mean(trueData{i,3}(3:4,idx),1);
    prediction_probability = trueData{i,5}; % row: early trials in tBlock, early trials in vBlock, late trials in tBlock, late trials in vBlock
    test_tBlock = mean(trueData{i,3}(1:2,:),1);
    test_vBlock = mean(trueData{i,3}(3:4,:),1);
    NumOfSessions = length(test_tBlock); 
    % Kendall's tau
    for b=1:B
        pred_early = prediction_probability(b,:);
        x_pred_early = repmat(2,length(pred_early),1);
        pred_late = prediction_probability(B+b,:);
        x_pred_late = repmat(3,length(pred_late),1);
        switch b
            case 1
                x_t = repmat(4,length(test_tBlock),1);
                x_v = repmat(1,length(test_tBlock),1);
%                 xtickNames = {'Light' 'Early' 'Late' 'Touch'};
                data_x = [x_v'; x_pred_early'; x_pred_late'; x_t'];
                data_y = [test_vBlock; pred_early; pred_late; test_tBlock];
                tail_type = 'right';
            case 2
                x_t = repmat(1,length(test_tBlock),1);
                x_v = repmat(4,length(test_tBlock),1);
%                 xtickNames = {'Touch' 'Early' 'Late' 'Light'};
                data_x = [x_t'; x_pred_early'; x_pred_late'; x_v'];
                data_y = [test_tBlock; pred_early; pred_late; test_vBlock];
                tail_type = 'left';   
        end
        subplot(2,B,(b-1)*2+1)
        scatter(x_pred_early,pred_early,'MarkerEdgeColor',[0.7 0.65 0.75],'Linewidth',1);hold on;
        scatter(x_pred_late,pred_late,'MarkerEdgeColor',[0.5 0.2 0.55],'Linewidth',1);
        scatter(x_t,test_tBlock,'MarkerEdgeColor',Colors_block{1});
        scatter(x_v,test_vBlock,'MarkerEdgeColor',Colors_block{2});
%         plot(data_y,'Color',[0.7 0.7 0.7])
        xlim([0.5 4.5]);ylim([0 1]);
        xticks(1:4);
        ylabel([Class1_name ' probability']);
        title(pred_names{b});
        xticklabels(xtickNames);
        set(gca, 'box','off', 'TickDir','out')
        for s=1:NumOfSessions
            y = data_y(:,s);
            x = data_x(:,s);
            idx = ~isnan(y); % remove NaN 
            [tau, pval] = corr(x(idx), y(idx), 'Type','Kendall', 'Tail', tail_type);
%             [tau, pval] = corr(x, y, 'Type','Kendall', 'Tail', tail_type);
            tau_ind_session(s) = tau;
            plot(x(idx),y(idx),'Color',[0.7 0.7 0.7])
        end
        tau_ind_avg = nanmean(tau_ind_session);
        switch b
            case 1
                text(3,0.1,['Tau = ' num2str(round(tau_ind_avg,2))]);
            case 2
                text(2,0.1,['Tau = ' num2str(round(tau_ind_avg,2))]);
        end
        clear tau_ind_session
    end
    
    %%% Bootstrapped prediction
    % Kendall's tau
    for boot=1:nBoot
%         idx = sum(isnan(bootstrpData{i,5}(:,:,boot)),1)==0; % Remove sessions without early or late trials
%         prediction_probability = bootstrpData{i,5}(:,idx,boot);
%         test_tBlock = mean(bootstrpData{i,3}(1:2,idx,boot),1);
%         test_vBlock = mean(bootstrpData{i,3}(3:4,idx,boot),1);
        prediction_probability = bootstrpData{i,5}(:,:,boot);
        test_tBlock = mean(bootstrpData{i,3}(1:2,:,boot),1);
        test_vBlock = mean(bootstrpData{i,3}(3:4,:,boot),1);
        NumOfSessions = length(test_tBlock); 
        % Kendall's tau
        for b=1:B
            pred_early =  prediction_probability(b,:);
            x_pred_early = repmat(2,length(pred_early),1);
            pred_late = prediction_probability(B+b,:);
            x_pred_late = repmat(3,length(pred_late),1);
            switch b
                case 1
                    x_t = repmat(4,length(test_tBlock),1);
                    x_v = repmat(1,length(test_tBlock),1);
                    data_x = [x_v'; x_pred_early'; x_pred_late'; x_t'];
                    data_y = [test_vBlock; pred_early; pred_late; test_tBlock];
                    tail_type = 'right';
                case 2
                    x_t = repmat(1,length(test_tBlock),1);
                    x_v = repmat(4,length(test_tBlock),1);
                    data_x = [x_t'; x_pred_early'; x_pred_late'; x_v'];
                    data_y = [test_tBlock; pred_early; pred_late; test_vBlock];
                    tail_type = 'left';   
            end
            for s=1:NumOfSessions
                y = data_y(:,s);
                x = data_x(:,s);
                idx = ~isnan(y); % remove NaN 
                [tau, pval] = corr(x(idx), y(idx), 'Type','Kendall', 'Tail', tail_type);
%                 [tau, pval] = corr(x, y, 'Type','Kendall', 'Tail', tail_type);
                tau_ind_session(s) = tau;
            end
            tau_ind_avg = nanmean(tau_ind_session);
            bootstrapped_tau_ind(b,boot) = tau_ind_avg;
            clear tau_ind_session
        end
    end
    
    % Plot
    test_tBlock = squeeze(mean(bootstrpData{i,3}(1:2,:,:),[1 2]));
    tBlock_sorted = sort(test_tBlock);
    tBlock_mean = mean(tBlock_sorted);
    tBlock_up95CI = tBlock_sorted(round(0.975*nBoot));
    tBlock_low95CI = tBlock_sorted(round(0.025*nBoot));
    test_vBlock = squeeze(mean(bootstrpData{i,3}(3:4,:,:),[1 2]));
    vBlock_sorted = sort(test_vBlock);
    vBlock_mean = mean(vBlock_sorted);
    vBlock_up95CI = vBlock_sorted(round(0.975*nBoot));
    vBlock_low95CI = vBlock_sorted(round(0.025*nBoot));
    for b = 1:B
        subplot(2,B,2*b)
        for g=1:2
            prediction = squeeze(nanmean(bootstrpData{i,5}((g-1)*B+b,:,:),[1,2]));
            prediction_sorted = sort(prediction);
            pred_mean = nanmean(prediction_sorted);
            pred_up95CI =  prediction_sorted(round(0.975*nBoot));
            pred_low95CI =  prediction_sorted(round(0.025*nBoot));

            errorbar(1+g,pred_mean,(pred_mean-pred_low95CI),(pred_up95CI-pred_mean),'o',...
                'Color', Colors_transition{g},'MarkerEdgeColor', Colors_transition{g},...
                'MarkerFaceColor',Colors_transition{g},'CapSize', 6, 'LineWidth', 0.5); hold on;
            switch b
                case 1
                    x_t = 4;v_t = 1;
%                     xtickNames = {'Light' 'Early' 'Late' 'Touch'};
                case 2
                    x_t = 1;v_t = 4;
%                     xtickNames = {'Touch' 'Early' 'Late' 'Light'};
            end
            if g ==2
                errorbar(x_t,tBlock_mean,(tBlock_mean-tBlock_low95CI),(tBlock_up95CI-tBlock_mean),'o',...
                    'Color', Colors_block{1},'MarkerEdgeColor', Colors_block{1},...
                    'MarkerFaceColor',Colors_block{1},'CapSize', 6, 'LineWidth', 0.5); hold on;

                errorbar(v_t,vBlock_mean,(vBlock_mean-vBlock_low95CI),(vBlock_up95CI-vBlock_mean),'o',...
                    'Color', Colors_block{2},'MarkerEdgeColor', Colors_block{2},...
                    'MarkerFaceColor',Colors_block{2},'CapSize', 6, 'LineWidth', 0.5); hold on;
            end
            
            xlim([0.5 4.5]);ylim([0 1]);
            xticks(1:4);
            xticklabels(xtickNames);
            ylabel([Class1_name ' probability']);
            title(pred_names{b})
            set(gca, 'box','off', 'TickDir','out')
        end
           
        % 95%CI tau            
        stat_sorted = sort(bootstrapped_tau_ind(b,:));
        stat_mean = round(mean(stat_sorted),3);
        stat_up95CI = round(stat_sorted(round(0.975*nBoot)),3);
        stat_low95CI = round(stat_sorted(round(0.025*nBoot)),3);
        text(0.8,0.1,['Tau = ' num2str(stat_mean) '(' num2str(stat_low95CI) ' ~ ' num2str(stat_up95CI) ')']);   
    end
    sgtitle(recSite_names{i})
%     % Save figure
%     predFigPath = fullfile(mainDir,['LDA_true_and_bootstrapped_prediction_', num2str(nCV),'CV_', num2str(nBoot),'nBoot_' recSite_names{i}]);
%     print(predFigPath,'-dpdf','-painters','-loose');
    clear bootstrapped_tau_ind
end

