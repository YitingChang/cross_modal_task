classdef Statistics
    % getShuffleData and getTestTrainTrials are modified dPCA functions.
    
    methods(Static)
      function [resampled_Xtrial] = bootstrp(Xtrial,numOfTrials,B,S,D)
            resampled_Xtrial = Xtrial;
            N = size(Xtrial,1);  
            %%%%%%%%% Bootstrap trials for each trial type %%%%%%%%
            for n = 1:N
                for b = 1:B
                    for s = 1:S
                        for d = 1:D
                            numOfTrials_n_tt = numOfTrials(n,b,s,d);
                            if numOfTrials_n_tt>0
                                Xtrial_n_tt = squeeze(Xtrial(n,b,s,d,:,1:numOfTrials_n_tt)); % remove NaN
                                isTrial = randi(numOfTrials_n_tt,numOfTrials_n_tt,1);
                                resampled_Xtrial(n,b,s,d,:,1:numOfTrials_n_tt) = Xtrial_n_tt(:,isTrial);
                            end
                        end
                    end
                end
            end   
        end
        
        function [resampled_Xtrial, resampled_numOfTrials] = Multi_bootstrp(Xtrial,numOfTrials,B,S,D)
            resampled_Xtrial = Xtrial;
            N = size(Xtrial,1);  
            %%%%%%%%% Bootstrap neurons %%%%%%%%
            isNeuron = randi(N,N,1);
            Xtrial_resampled_n = Xtrial(isNeuron,:,:,:,:,:); 
            numOfTrials_resampled_n = numOfTrials(isNeuron,:,:,:); 

            %%%%%%%%% Bootstrap trials for each trial type %%%%%%%%
            for n = 1:N
                for b = 1:B
                    for s = 1:S
                        for d = 1:D
                            numOfTrials_n_tt = numOfTrials_resampled_n(n,b,s,d);
                            if numOfTrials_n_tt>0
                                Xtrial_n_tt = squeeze(Xtrial_resampled_n(n,b,s,d,:,1:numOfTrials_n_tt)); % remove NaN
                                isTrial = randi(numOfTrials_n_tt,numOfTrials_n_tt,1);
                                resampled_Xtrial(n,b,s,d,:,1:numOfTrials_n_tt) = Xtrial_n_tt(:,isTrial);
                                resampled_numOfTrials(n,b,s,d) = numOfTrials_n_tt;
                            end
                        end
                    end
                end
            end   
        end
        
        function [resampled_Xtrial, resampled_numOfTrials] = Multi_bootstrp_correct(Xtrial,numOfTrials,S,C)
                resampled_Xtrial = Xtrial;
                N = size(Xtrial,1);  
                %%%%%%%%% Bootstrap neurons %%%%%%%%
                isNeuron = randi(N,N,1);
                FR_resampled_n = Xtrial(isNeuron,:,:,:,:); 
                numOfTrials_resampled_n = numOfTrials(isNeuron,:,:); 

                %%%%%%%%% Bootstrap trials for each trial type %%%%%%%%
                for n = 1:N
                    for s = 1:S
                        for c = 1:C
                            numOfTrials_n_tt = numOfTrials_resampled_n(n,s,c);
                            FR_n_tt = squeeze(FR_resampled_n(n,s,c,:,1:numOfTrials_n_tt)); % remove NaN
                            isTrial = randi(numOfTrials_n_tt,numOfTrials_n_tt,1);
                            resampled_Xtrial(n,s,c,:,1:numOfTrials_n_tt) = FR_n_tt(:,isTrial);
                            resampled_numOfTrials(n,s,c) = numOfTrials_n_tt;
                        end
                    end
                end   
        end
        
        function firingRates = getDataArray4dPCA(firingRates, recSite, timeBin, B)
            % timeBin = [start, end]
            block_decision = {[1,3]; [2,3]};
            recSite_idx = find(strcmp(firingRates(:,1),recSite),1);
            Xtrial = firingRates{recSite_idx,2}(:,:,:,:,timeBin(1):timeBin(2),:);
            X = firingRates{strcmp(firingRates(:,1),recSite),3}(:,:,:,:,timeBin(1):timeBin(2),:);
            trialNum = firingRates{strcmp(firingRates(:,1),recSite),4};
            % normalization
            X_mean = nanmean(X(:,:),2);
            X_range = range(X(:,:),2);
            normalized_Xtrial = (Xtrial - X_mean)./(X_range + 0.01);  
            normalized_X = (X - X_mean)./(X_range + 0.01);  
            for b = 1:B
                firingRates{recSite_idx,5+ (b-1)*3} = squeeze(normalized_Xtrial(:,b,:,block_decision{b},:,:));
                firingRates{recSite_idx,6+ (b-1)*3} = squeeze(normalized_X(:,b,:,block_decision{b},:));
                firingRates{recSite_idx,7+ (b-1)*3} = squeeze(trialNum(:,b,:,block_decision{b}));
            end
        end            
        
%         function [XShuffle, XtrialShuffle] = getShuffledData(Xtrial,numOfTrials)
%             % shuffling PSTHs inside each neuron (preserving time and numOfTrials)
%             
%             D = size(numOfTrials,1);
%             dim = size(Xtrial);
%             T = dim(end-1);
%             maxTrialN = size(Xtrial, ndims(Xtrial));
%             numCond = prod(dim(2:end-2));
% 
%             % find missing trials
%             numOfTrialsCond = reshape(numOfTrials, D, []);
%             trialsMissing = zeros(D, maxTrialN*numCond);
%             for n=1:D
%                 this = zeros(numCond, maxTrialN);
%                 for c=1:numCond
%                     this(c,:) = [zeros(1,numOfTrialsCond(n,c)) ones(1,maxTrialN-numOfTrialsCond(n,c))];
%                 end
%                 trialsMissing(n,:) = this(:);
%             end
% 
%             % collapsing conditions
%             orderDim = 1:ndims(Xtrial);
%             orderDim(end-1:end) = orderDim([end end-1]);
%             XtrialCond = permute(Xtrial, orderDim); % time shifted to the end
%             XtrialCond = reshape(XtrialCond, D, [], T);
% 
%             % shuffling PSTHs inside each neuron (preserving time and numOfTrials)
%             XtrialCondShuffle = zeros(size(XtrialCond));
% 
%             for n = 1:D
%                 presentTrials = find(trialsMissing(n,:) == 0);
%                 shuffledOrder = presentTrials(randperm(length(presentTrials)));
%                 XtrialCondShuffle(n,presentTrials,:) = XtrialCond(n,shuffledOrder,:);
%             end
% 
%             XtrialShuffle = permute(reshape(XtrialCondShuffle, dim(orderDim)), ...
%                 orderDim);
%             clear XtrialCondShuffle
% 
%             XShuffle = sum(XtrialShuffle, ndims(XtrialShuffle));
%             XShuffle = bsxfun(@times, XShuffle, 1./numOfTrials);
%         end
        
        function [XShuffle, XtrialShuffle] = getShuffledData(Xtrial,numOfTrials)
            % shuffling PSTHs inside each neuron (preserving time and numOfTrials)
            % Xtrial: NxSxCxTxK
            dim = size(Xtrial);
            N = dim(1); % number of neurons
            T = dim(end-1); % number of time bins
            maxTrialNum = dim(end); % maximum trial numbers
            numCond = prod(dim(2:end-2)); % number of conditions

            % find missing trials
            numOfTrialsCond = reshape(numOfTrials, N, []);
            trialsMissing = zeros(N, maxTrialNum*numCond);
            for n=1:N
                this = zeros(numCond, maxTrialNum);
                for c=1:numCond
                    this(c,:) = [zeros(1,numOfTrialsCond(n,c)) ones(1,maxTrialNum-numOfTrialsCond(n,c))];
                end
                trialsMissing(n,:) = this(:);
            end

            % collapsing conditions
            orderDim = 1:ndims(Xtrial);
            orderDim(end-1:end) = orderDim([end end-1]);
            XtrialCond = permute(Xtrial, orderDim); % time shifted to the end
            XtrialCond = reshape(XtrialCond, N, [], T);

            % shuffling PSTHs inside each neuron (preserving time and numOfTrials)
            XtrialCondShuffle = zeros(size(XtrialCond));

            for n = 1:N
                presentTrials = find(trialsMissing(n,:) == 0);
                shuffledOrder = presentTrials(randperm(length(presentTrials)));
                XtrialCondShuffle(n,presentTrials,:) = XtrialCond(n,shuffledOrder,:);
            end

            XtrialShuffle = permute(reshape(XtrialCondShuffle, dim(orderDim)), ...
            orderDim);
            clear XtrialCondShuffle

            XShuffle = sum(XtrialShuffle, ndims(XtrialShuffle));
            XShuffle = bsxfun(@times, XShuffle, 1./numOfTrials);
        end
        
        
        
        function [Xtest, Xtrainfull] = getTestTrainTrials(Xtrial, numOfTrials, varargin)
             
            % For each trial condition, one trial is selected as a test trial and the
            % rest of trials are training trials.

            dim = size(Xtrial);

            neuronsConditions = numOfTrials(:);
            testTrials = ceil(rand([length(neuronsConditions) 1]) .* neuronsConditions);

            ind = reshape(testTrials, size(numOfTrials));
            ind = bsxfun(@times, ones(dim(1:end-1)), ind);
            ind = ind(:);

            indtest = sub2ind([prod(dim(1:end-1)) dim(end)], (1:prod(dim(1:end-1)))', ind);

            Xtest = Xtrial(indtest);
            Xtest = reshape(Xtest, dim(1:end-1));

            if nargout > 1
                Xtrainfull = Xtrial;
                Xtrainfull(indtest) = nan;
            end
        end
        
        function [Ztest, Ztrain] = getTestTrainTrials_CV(Xtrial, X, numOfTrials,S,D,C,nCV,varargin)
            % inputs: 
            % Xtrial: firing rates per trial
            % X     : mean firing rates
            % numOfTrials
            
            % For each trial condition, 10% of trials are selected as test trials and the
            % rest of trials are training trials.
            
            % Normalization
            X_mean = mean(X(:,:),2);
            X_range = range(X(:,:),2);
            Ztrial = (Xtrial - X_mean)./(X_range);   
            
            N = size(X,1);
            T = size(X,length(size(X))); % time is the last dimension
            for n = 1:N
                for s = 1:S
                    for d = 1:D
                        for c = 1:C
                            numOfTrials_tt = numOfTrials(n,s,d,c);
                            X_n_tt = squeeze(Ztrial(n,s,d,c,:,1:numOfTrials_tt)); % remove NaN
                            isTest = randperm(numOfTrials_tt) <= round(numOfTrials_tt/nCV,0);
                            Ztest(n,s,d,c,1:T)= mean(X_n_tt(:,isTest),2);
                            Ztrain(n,s,d,c,1:T)= mean(X_n_tt(:,~isTest),2);
                        end
                    end
                end
            end   
        end
        
        function [Ztest, Ztrain] = getTestTrainTrials_CV_correct(Xtrial, X, numOfTrials,S,C,nCV,varargin)
            % inputs: 
            % Xtrial: firing rates per trial
            % X     : mean firing rates
            % numOfTrials
            
            % For each trial condition, 10% of trials are selected as test trials and the
            % rest of trials are training trials.
            
            % Normalization
            X_mean = mean(X(:,:),2);
            X_range = range(X(:,:),2);
            Ztrial = (Xtrial - X_mean)./(X_range);   
            
            N = size(X,1);
            T = size(X,length(size(X))); % time is the last dimension
            for n = 1:N
                for s = 1:S
                    for c = 1:C
                        numOfTrials_tt = numOfTrials(n,s,c);
                        X_n_tt = squeeze(Ztrial(n,s,c,:,1:numOfTrials_tt)); % remove NaN
                        isTest = randperm(numOfTrials_tt) <= round(numOfTrials_tt/nCV,0);
                        Ztest(n,s,c,1:T)= mean(X_n_tt(:,isTest),2);
                        Ztrain(n,s,c,1:T)= mean(X_n_tt(:,~isTest),2);
                    end
                end
            end   
        end
    end
end