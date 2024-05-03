%% Lick probability for individual mice
% Select lick responses based on brain regions of interest
load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav')
Genotypes = {'PV-Cre-Ai32'}; 
% Genotypes = {'PV-Cre-Ai32' 'VGAT-ChR2' 'all'};
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
% inhSites = {'ALM'}; 

save_path = 'E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\Figures_lick_prob\';

B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 

nboot = 10000; % number for bootstrapping
colors_lick = {'b' 'r'};
titles_tt = {'TT' 'TV'; 'VT' 'VV'};
titles_catch = {'T\_catch' 'V\_catch'};

for g = 1:length(Genotypes)
    Genotype = Genotypes{g};
    for e = 1:length(inhSites) 
        inhSite = inhSites{e};
        if g<3
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite)& strcmp(photoInh_all.Genotype, Genotype), :);
        else
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite), :); % pooled genotypes
        end
     
        mouseNames = unique(u.MouseName);

        for mouse = 1:length(mouseNames)
            mouseName = mouseNames{mouse};
            m = u(strcmp(u.MouseName, mouseName), :);
            % Get lick responses for trials with stimuli
            % Bootstrap 
            for b=1:B
                for s=1:2
                    for l=1:L
                        for o=1:O
                            v=[];
                            % Concatenation
                            for i = 1:height(m)
                                new_v = m{i,5}{1,1}{b,s,l,o}; 
                                v = [v; new_v];
                            end
                            if length(v) >1
                                [ci bootstat]= bootci(nboot,@mean,v);  
                                MEAN = mean(bootstat);
                                data_bootstrapped(b,s,l,o,1) = MEAN;
                                data_bootstrapped(b,s,l,o,2) = ci(1);
                                data_bootstrapped(b,s,l,o,3) = ci(2);
                            else
                                data_bootstrapped(b,s,l,o,1) = v;
                                data_bootstrapped(b,s,l,o,2) = 0;
                                data_bootstrapped(b,s,l,o,3) = 0;
                            end
                        end
                    end
                end
            end
            % Get lick responses for trials without stimuli
            for b=1:B
                for s=3:4
                    for l=1:L
                        for o= 1:O
                            v=[];
                            % Concatenation
                            for i = 1:height(m)
                                new_v = m{i,5}{1,1}{b,s,l,o}; 
                                v = [v; new_v];
                            end
                            if length(v)>1
                                [ci bootstat]= bootci(nboot,@mean,v);  
                                MEAN = mean(bootstat);
                                data_bootstrapped(b,s,l,o,1) = MEAN;
                                data_bootstrapped(b,s,l,o,2) = ci(1);
                                data_bootstrapped(b,s,l,o,3) = ci(2);
                            elseif length(v)== 1
                                data_bootstrapped(b,s,l,o,1) = v;
                                data_bootstrapped(b,s,l,o,2) = 0;
                                data_bootstrapped(b,s,l,o,3) = 0;
                            else
                                data_bootstrapped(b,s,l,o,1) = NaN;
                                data_bootstrapped(b,s,l,o,2) = NaN;
                                data_bootstrapped(b,s,l,o,3) = NaN;                             
                            end
                        end
                    end
                end
            end
            % plot for individual mouse
            figure('Position', [0,0, 600, 600]); % set figure size
            for b=1:B
                for s=1:2    
                    for l=1:2
                        for o=1:O
                            subplot(2,3,(b-1)*3+s+1)
                            x = o;
                            y1 = data_bootstrapped(b,s,l,o,1);
                            y2 = data_bootstrapped(b,s,l,o,2);
                            y3 = data_bootstrapped(b,s,l,o,3);
                            errorbar(x,y1,(y1-y2),(y3-y1),'o','Color',colors_lick{l},...
                            'MarkerSize',6,'MarkerFaceColor',colors_lick{l}, 'MarkerEdgeColor', colors_lick{l})
                            hold on
                        end
                    end
                    set(gca, 'box','off','TickDir','out');
                    ylim([0 1]);
                    ax = gca; ax.YAxis.Visible = 'off';
                    xlim([0.5 3.5]);  
                    xtickNames = {'no manipulation' 'pre-inhibition' 'post-inhibition'};
                    xticklabels(xtickNames);
                    xtickangle(45)
                    title(titles_tt(b,s))
                end
            end
            % catch trial
            for b=1:B
                for s=3:4
                    for l=1:2
                        if s==3
                            o_idx = [1,2]; x_val = {1 2 []};
                        else
                            o_idx = [1,3]; x_val = {3 [] 4};
                        end
                        for o= o_idx
                            subplot(2,3,(b-1)*3+1)
                            x = x_val{o};
                            y1 = data_bootstrapped(b,s,l,o,1);
                            y2 = data_bootstrapped(b,s,l,o,2);
                            y3 = data_bootstrapped(b,s,l,o,3);
                            errorbar(x,y1,(y1-y2),(y3-y1),'o','Color',colors_lick{l},...
                            'MarkerSize',6,'MarkerFaceColor',colors_lick{l}, 'MarkerEdgeColor', colors_lick{l})
                            hold on
                        end
                    end
                    set(gca, 'box','off','TickDir','out');
                    ylim([0 1]);
                    ylabel('lick probability');
                    xlim([0.5 4.5]); 
                    xticks(1:4)
                    xticklabels({'short\_ITI' 'long\_ITI' 'short\_laser only' 'long\_laser only'})
                    xtickangle(45)
                end
                title(titles_catch{b})
            end
            sgtitle([mouseName '\_' inhSite], 'FontSize',20);
            % Save fig
            fig = gcf;
            inhFigPath = fullfile(save_path,[mouseName,'_' inhSite '.pdf']);
            print(inhFigPath,'-dpdf','-painters','-loose');
            clear data_bootstrapped
        end
    end
end

%% compare to no manipulation
% Multilevel_bootstrap
% resample mice with replacement
% resample session with replacement
% resample trials with replacement 

load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 

Genotypes = {'PV-Cre-Ai32'}; 
B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 
% catch trial: all (s=3), short (s=4), long (s=5)
% o=1: ITI of catch trials; o=3: laser catch trials

nboot = 10000; % number for bootstrapping

for g = 1:length(Genotypes)
    Genotype = Genotypes{g};
    for e = 1:length(inhSites) 
        inhSite = inhSites{e};
        if g<3
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite)& strcmp(photoInh_all.Genotype, Genotype), :);
        else
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite), :); % pooled genotypes
        end
        mouseNames = unique(u.MouseName);
        rng(7);
%             rng(test(t)); % control random number generation
        for n = 1:nboot
            resampled_mice = randsample(mouseNames, length(mouseNames), true); 
            for mouse = 1:length(resampled_mice)
                m = u(strcmp(u.MouseName, resampled_mice{mouse}), :);
                % resample session
                for session = 1:height(m)
                    v = m{randi(height(m)),5}{1,1};   
                    % Trials with stimuli
                    for b=1:B
                        for s=1:2
                            for l=1:L-1
                                for o=1:O
                                    % resample trials
                                    lick_tt = v{b,s,l,o};
                                    resampled_lick_tt = randsample(lick_tt, length(lick_tt), true); 
                                    lick_prob(b,s,l,o) = mean(resampled_lick_tt); % lick probability of each trial type
                                    if o > 1
                                        delta_lick_prob_session(session,b,s,l,o-1) = lick_prob(b,s,l,o) - lick_prob(b,s,l,1);
                                        % difference in lick probability (inhibition - no inhibition)
                                    end
                                end
                            end
                        end
                    end
                    % Trials without stimuli
                    for b=1:B
                        for s=3:4
                            if s==3
                                o_idx = [1 2];
                            else
                                o_idx = [1 3];
                            end
                            for l=1:L-1
                                for o= o_idx
                                    % resample trials
                                    lick_tt = v{b,s,l,o};
                                    resampled_lick_tt = randsample(lick_tt, length(lick_tt), true); 
                                    lick_prob(b,s,l,o) = mean(resampled_lick_tt); % lick probability of each trial type
                                    if o > 1
                                        delta_lick_prob_session(session,b,s,l,1) = lick_prob(b,s,l,o) - lick_prob(b,s,l,1);
                                        % difference in lick probability (inhibition - no inhibition)
                                    end
                                end
                            end
                        end
                    end
                    clear lick_prob
                end
                delta_lick_prob_mouse(mouse,1:B,1:S,1:(L-1),1:(O-1)) = mean(delta_lick_prob_session,1); % average across sessions
                clear delta_lick_prob_session
            end
            delta_lick_prob(n,1:B,1:S,1:(L-1),1:(O-1)) = mean(delta_lick_prob_mouse,1); % average across mice
            clear delta_lick_prob_mouse
        end

        photoInh_bootstrapped{(g-1)*length(inhSites)+e,1} = Genotype;
        photoInh_bootstrapped{(g-1)*length(inhSites)+e,2} = inhSite;
        photoInh_bootstrapped{(g-1)*length(inhSites)+e,3} = delta_lick_prob;
        clear delta_lick_prob
    end
end
% % Convert to a table
% VarNames = {'Genotype' 'inhSite' 'bootstrapping'};
% photoInh_bootstrap = cell2table(photoInh_bootstrap,'VariableNames', VarNames); 

% Save 
photoInhPaths = fullfile('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control', ... 
    'photoinhibition_bootstrapped_PV-Cre-Ai32');

save(photoInhPaths, 'photoInh_bootstrapped'); 
% clear photoInh_bootstrapped 

%% Plotting
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
Genotypes = {'PV-Cre-Ai32'}; 
sgtitles_tt = {'TT' 'TV'; 'VT' 'VV'};
sgtitles_catch = {'short catch' 'long catch'};
titles_tt = {'pre-inhibition' 'post-inhibition'};
titles_block = {'resond-to-touch' 'respond-to-light'};
ylabels = {'delta right lick prob' 'delta left lick prob'};
colors_tt = {[0 0 1] [1 0 1]; [0 1 1] [1 0 0]}; 
% colors_block = {[0.25 0.5 0.05] [1 0.4 0.16]};
colors_block = {'b' 'b'};
xtickNames = {'sham' 'S1/S2' 'MM' 'ALM' 'AMM' 'Prt'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\Figures_delta_lick_prob';
 
B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 
% catch trial: all (s=3), short (s=4), long (s=5)
% o=1: ITI of catch trials; o=3: laser catch trials

nboot = 10000; % number for bootstrapping
% for one gynotype 

% Trials with stimuli
for b=1:B
    for s=1:2   
        figure('Position', [0,0, 600, 600]); % set figure size
        for l=1:L-1
            for o=1:O-1
                subplot(2,2,(l-1)*(L-1)+o)
                for e = 1:length(inhSites) 
                    delta_lick_prob_bootstrapped_tt = photoInh_bootstrapped{e,3}(:,b,s,l,o);
                    sorted_delta_lick_prob = sort(delta_lick_prob_bootstrapped_tt);   
                    y1 = mean(sorted_delta_lick_prob); 
                    y2 = sorted_delta_lick_prob(nboot*0.025); 
                    y3 = sorted_delta_lick_prob(nboot*0.975); 
                    x = e;
                    errorbar(x,y1,(y1-y2),(y3-y1),'o','Color',colors_tt{b,s},...
                    'MarkerSize',6,'MarkerFaceColor',colors_tt{b,s}, 'MarkerEdgeColor',colors_tt{b,s})
                    hold on
                    CI(o,e,1) = y2;
                    CI(o,e,2) = y3;
                end
                set(gca, 'box','off','TickDir','out');
                ylim([-1 0.6]);
                ylabel(ylabels{l});
                yline(0, '--k');
                xlim([0.5 (length(inhSites)+0.5)]);  
                xticks([1:length(inhSites)]);
                xticklabels(xtickNames);
                xtickangle(45)
                title(titles_tt{o})
            end
        end
        sgtitle(sgtitles_tt(b,s))
%         save fig       
%         fig = gcf;
%         inhFigPath = fullfile(save_path,['CM_inhibition_' sgtitles_tt{b,s} '.pdf']);
%         print(inhFigPath,'-dpdf','-painters','-loose');   
    end
end

% Trials without stimuli

for s=3:4
    figure('Position', [0,0, 600, 600]); % set figure size
    for b=1:B 
        for l=1:2   
            subplot(2,2,(l-1)*(L-1)+b)
            for e = 1:length(inhSites) 
                isInhSites = cellfun(@(x) strcmp(x,inhSites{e}),photoInh_bootstrapped(:,2));
                delta_lick_prob_bootstrapped_tt = photoInh_bootstrapped{isInhSites,3}(:,b,s,l,1);
                sorted_delta_lick_prob = sort(delta_lick_prob_bootstrapped_tt);   
                y1 = mean(sorted_delta_lick_prob); 
                y2 = sorted_delta_lick_prob(nboot*0.025); 
                y3 = sorted_delta_lick_prob(nboot*0.975); 
                x = e;
                errorbar(x,y1,(y1-y2),(y3-y1),'o','Color',colors_block{b},...
                'MarkerSize',6,'MarkerFaceColor',colors_block{b}, 'MarkerEdgeColor',colors_block{b})
                hold on
            end
            set(gca, 'box','off','TickDir','out');
            ylim([-1 0.5]);
            set(gca,'Ytick',-1:0.5:0.5); 
            ylabel(ylabels{l});
            yline(0, '--k');
            xlim([0.5 (length(inhSites)+0.5)]);  
            xticks([1:length(inhSites)]);
            xticklabels(xtickNames);
            xtickangle(45)
            title(titles_block{b})
        end
    end
    sgtitle(sgtitles_catch{s-2})
    % save fig       
    fig = gcf;
    inhFigPath = fullfile(save_path,['CM_inhibition_' sgtitles_catch{s-2} '.pdf']);
    print(inhFigPath,'-dpdf','-painters','-loose'); 
end

%% sensitivity index
% respond-to-touch block: sensitivity = tHit - vFA_right_lick
% respond-to-light block: sensitivity = vHit - tFA_left_lick

% Select lick responses based on brain regions of interest
load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav')
Genotypes = {'PV-Cre-Ai32'}; 
% Genotypes = {'PV-Cre-Ai32' 'VGAT-ChR2' 'all'};
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
% inhSites = {'ALM'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\Figures_sensitivity\';

B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 
% catch trial: all (s=3), short (s=4), long (s=5)
% o=1: ITI of catch trials; o=3: laser catch trials
tt = {[1 1 1] [1 2 1]; [2 2 2] [2 1 2]}; % tHit vFA; vHit tFA
% colors = {[0 0 0] [0.6 0.9 0.5] [0.25 0.5 0.05]}; % no, pre- ,post- inhibition
blockNames = {'respond-to-touch' 'respond-to-light'};
optoNames = {'control' 'pre-inh.' 'post-inh'};
sgtitles = {'sham' 'S1\_S2' 'MM' 'ALM' 'AMM' 'Prt'};
ylabels = {'Hit rate' 'FA rate' 'sensitivity'};
optoColors = {'k' 'b' 'r'};
for g = 1:length(Genotypes)
    Genotype = Genotypes{g};
    for e = 1:length(inhSites) 
        inhSite = inhSites{e};
        if g<3
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite)& strcmp(photoInh_all.Genotype, Genotype), :);
        else
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite), :); % pooled genotypes
        end
        mouseNames = unique(u.MouseName);
        for mouse = 1:length(mouseNames)
            mouseName = mouseNames{mouse};
            m = u(strcmp(u.MouseName, mouseName), :);
            for session = 1:height(m)
                v = m{session,5}{1,1};
                for b=1:B
                    for o=1:O
                        Hit = v{tt{b,1}(1),tt{b,1}(2),tt{b,1}(3),o}; 
                        pHit = mean(Hit);% hit rate
                        FA = v{tt{b,2}(1),tt{b,2}(2),tt{b,2}(3),o}; 
                        pFA = mean(FA);% false alarm rate
                        sensitivity = pHit - pFA;
%                         dpri = dprime(pHit,pFA,length(Hit),length(FA));
                        data_session(1,session,b,o) = pHit;
                        data_session(2,session,b,o) = pFA;
                        data_session(3,session,b,o) = sensitivity;
                    end
                end
            end
            data_mouse(1:3,mouse,1:B,1:O) = mean(data_session,2);
            clear data_session
        end
        % plot all mice
        figure('Position', [0,0, 600, 900]); % set figure size
        for k=1:3
            for b=1:B
                subplot(3,B,(k-1)*B+b)
                for o=1:O
                    y = data_mouse(k,:,b,o);
                    x = repmat(o,length(y),1);
                    scatter(x,y,'b');hold on;
                end
                line = squeeze(data_mouse(k,:,b,:))';
                plot(line,'Color',[0.7 0.7 0.7])  
                title(blockNames{b})
                xlim([0.5 3.5])
                xticklabels(optoNames)
                ylabel(ylabels{k})
                if k == 3
                    ylim([-1 1])
                else
                    ylim([0 1])
                end
            end
        end
        sgtitle(sgtitles{e})
        % Save fig
        fig = gcf;
        inhFigPath = fullfile(save_path,['sensitivity_' inhSite '.pdf']);
        print(inhFigPath,'-dpdf','-painters','-loose');
        clear dp_mouse pHit_mouse pFA_mouse
        
        figure('Position', [0,0, 600, 400]); % set figure size
        for b=1:B
            for o=1:O
                subplot(2,O,(b-1)*O+o)
                for k = 1:2
                    y = data_mouse(k,:,b,o);
                    x = repmat(k,length(y),1);
                    scatter(x,y,'b');hold on;
                end
                line = squeeze(data_mouse(1:2,:,b,o));
                plot(line,'Color',[0.7 0.7 0.7])
                set(gca, 'box','off','TickDir','out');
                ylim([0 1]);
                xlim([0.5 2.5])
                xticks(1:2); 
                xticklabels({'Hit' 'FA'})
                title(optoNames{o})
                if o == 1
                    ylabel({blockNames{b};'p(lick)'})
                else                
                    ax = gca; ax.YAxis.Visible = 'off';
                end  
            end
        end
        sgtitle(sgtitles{e}) 
        % Save fig
        fig = gcf;
        inhFigPath = fullfile(save_path,['Hit_FA_' inhSite '.pdf']);
        print(inhFigPath,'-dpdf','-painters','-loose');
        
        clear data_mouse
    end
end

%% Change in sensitivity
% no manipulation vs. pre- or post- inhibition
% Multilevel_bootstrap
% resample mice with replacement
% resample session with replacement
% resample trials with replacement 
load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav')
Genotypes = {'PV-Cre-Ai32'}; 
% Genotypes = {'PV-Cre-Ai32' 'VGAT-ChR2' 'all'};
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
% inhSites = {'ALM'}; 

B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 
tt = {[1 1 1] [1 2 1]; [2 2 2] [2 1 2]}; % tHit vFA; vHit tFA
nboot = 10000; % number for bootstrapping     

for g = 1:length(Genotypes)
    Genotype = Genotypes{g};
    for e = 1:length(inhSites) 
        inhSite = inhSites{e};
        if g<3
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite)& strcmp(photoInh_all.Genotype, Genotype), :);
        else
            u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite), :); % pooled genotypes
        end
        mouseNames = unique(u.MouseName);
        rng(4);
        for n = 1:nboot
            resampled_mice = randsample(mouseNames, length(mouseNames), true); 
            for mouse = 1:length(resampled_mice)
                m = u(strcmp(u.MouseName, resampled_mice{mouse}), :);
                % resample session
                for session = 1:height(m)
                    v = m{randi(height(m)),5}{1,1};   
                    for b=1:B
                        for o=1:O
                            % resample hit trials      
                            Hit = v{tt{b,1}(1),tt{b,1}(2),tt{b,1}(3),o}; 
                            resampled_Hit = randsample(Hit, length(Hit), true); 
                            pHit = mean(resampled_Hit);% hit rate
                            % resample FA trials
                            FA = v{tt{b,2}(1),tt{b,2}(2),tt{b,2}(3),o}; 
                            resampled_FA = randsample(FA, length(FA), true);
                            pFA = mean(FA);% false alarm rate
                            sensitivity = pHit - pFA;
%                             dpri = dprime(pHit,pFA,length(Hit),length(FA));
                            st_session(session,b,o) = sensitivity; 
                            if o > 1
                                delta_st_session(session,b,o-1) = sensitivity - st_session(session,b,1);
                            end
                        end
                    end
                end
                delta_st_mouse(mouse,1:B,1:O-1) = mean(delta_st_session,1); 
                clear st_session delta_st_session
            end
            delta_sensitivity(n,1:B,1:O-1) = mean(delta_st_mouse,1);  
            clear delta_st_mouse
        end
        delta_sensitivity_bootstrapped{(g-1)*length(inhSites)+e,1} = Genotype;
        delta_sensitivity_bootstrapped{(g-1)*length(inhSites)+e,2} = inhSite;
        delta_sensitivity_bootstrapped{(g-1)*length(inhSites)+e,3} = delta_sensitivity;
        clear delta_sensitivity
    end
end
% Save 
photoInhPaths = fullfile('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control', ... 
    'delta_sensitivity_bootstrapped_PV-Cre-Ai32');

save(photoInhPaths, 'delta_sensitivity_bootstrapped'); 
        


%% Plot
% analysis setting
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
Genotypes = {'PV-Cre-Ai32'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\Figures_delta_sensitivity\';
B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 
nboot = 10000; % number for bootstrapping  

% figure setting
ylabels = {{'pre-inhibition'; 'delta sensitivity'} {'post-inhibition'; 'delta sensitivity'}};
titles_block = {'resond-to-touch' 'respond-to-light'};
xtickNames = {'sham' 'S1/S2' 'MM' 'ALM' 'AMM' 'Prt'};
colors_inhibition = {[0.65 0.9 0.5] [0.25 0.5 0.05]};

figure('Position', [0,0, 600, 600]); % set figure size
for b=1:B
    for o=1:O-1
        subplot(2,2,(o-1)*B+b)
        for e = 1:length(inhSites) 
            delta_sensitivity_bootstrapped_tt = delta_sensitivity_bootstrapped{e,3}(:,b,o);
            sorted_delta_sensitivity = sort(delta_sensitivity_bootstrapped_tt);   
            y1 = mean(sorted_delta_sensitivity); 
            y2 = sorted_delta_sensitivity(nboot*0.025); 
            y3 = sorted_delta_sensitivity(nboot*0.975); 
            x = e;
            errorbar(x,y1,(y1-y2),(y3-y1),'o','Color',colors_inhibition{o},...
            'MarkerSize',6,'MarkerFaceColor',colors_inhibition{o}, 'MarkerEdgeColor',colors_inhibition{o})
            hold on
        end
        set(gca, 'box','off','TickDir','out');
        ylim([-1 0.5]);
%         ylim([-3.5 1]);
        ylabel(ylabels{o});
        yline(0, '--k');
        xlim([0.5 (length(inhSites)+0.5)]);  
        xticks([1:length(inhSites)]);
        xticklabels(xtickNames);
        xtickangle(45)
        title(titles_block{b})
    end
end
% save fig       
fig = gcf;
inhFigPath = fullfile(save_path,'CM_inhibition_delta_sensitivity.pdf');
print(inhFigPath,'-dpdf','-painters','-loose'); 

%% Figure 6b-d
% delta sensitivity in the respond-to-touch blocks
% analysis setting
load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\delta_sensitivity_bootstrapped_PV-Cre-Ai32')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\Figures_delta_sensitivity\';
B=2; % respond-to-touch and respond-to-light blocks
S=4; % tactile, visual, short laser catch, long laser catch
L=3; % right, left, no lick
O=3; % no, pre-, post- inhibition 
nboot = 10000; % number for bootstrapping

% figure setting
figure('Position', [0,0, 700, 300]); % set figure size
xtickNames = {'sham' 'S1/S2' 'MM' 'ALM' 'AMM' 'Prt'};
titles = {'pre-inhibition' 'post-inhibition'};

for o=1:O-1
    subplot(1,2,o)
    for e = 1:length(inhSites) 
        delta_sensitivity_bootstrapped_tt = delta_sensitivity_bootstrapped{e,3}(:,2,o);
        sorted_delta_sensitivity = sort(delta_sensitivity_bootstrapped_tt);   
        y1 = mean(sorted_delta_sensitivity); 
        y2 = sorted_delta_sensitivity(nboot*0.025); 
        y3 = sorted_delta_sensitivity(nboot*0.975); 
        x = e;
        errorbar(x,y1,(y1-y2),(y3-y1),'o','Color','g',...
        'MarkerSize',6,'MarkerFaceColor','g', 'MarkerEdgeColor','g')
        hold on
        CI(o,e,1) = y2;
        CI(o,e,2) = y3;
    end
    set(gca, 'box','off','TickDir','out');
    ylim([-1 0.5]);
    ylabel('delta sensitivity');
    yline(0, '--k');
    xlim([0.5 (length(inhSites)+0.5)]);  
    xticks([1:length(inhSites)]);
    xticklabels(xtickNames);
    xtickangle(45)
    title(titles{o})
end
% save fig       
fig = gcf;
inhFigPath = fullfile(save_path,'CM_inhibition_delta_sensitivity_vBlock.pdf');
print(inhFigPath,'-dpdf','-painters','-loose'); 

%% Figure 6i
% delta right lick probability of short laser catch trials in respond-to-touch blocks
% analysis setting
load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoinhibition_bootstrapped_PV-Cre-Ai32')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'};
save_path = 'E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\Figures_delta_lick_prob';
nboot = 10000; % number for bootstrapping

% figure setting
figure('Position', [0,0, 300, 300]); % set figure size
xtickNames = {'sham' 'S1/S2' 'MM' 'ALM' 'AMM' 'Prt'}; 
for e = 1:length(inhSites) 
    isInhSites = cellfun(@(x) strcmp(x,inhSites{e}),photoInh_bootstrapped(:,2));
    delta_lick_prob_bootstrapped_tt = photoInh_bootstrapped{isInhSites,3}(:,1,4,1,1);
    sorted_delta_lick_prob = sort(delta_lick_prob_bootstrapped_tt);   
    y1 = mean(sorted_delta_lick_prob); 
    y2 = sorted_delta_lick_prob(nboot*0.025); 
    y3 = sorted_delta_lick_prob(nboot*0.975); 
    x = e;
    errorbar(x,y1,(y1-y2),(y3-y1),'o','Color','b',...
    'MarkerSize',6,'MarkerFaceColor','b', 'MarkerEdgeColor','b')
    hold on
    CI(e,1) = y2;
    CI(e,2) = y3;
end
set(gca, 'box','off','TickDir','out');
ylim([-1 0.5]);
set(gca,'Ytick',-1:0.5:0.5); 
ylabel('delta probability of right lick');
yline(0, '--k');
xlim([0.5 (length(inhSites)+0.5)]);  
xticks([1:length(inhSites)]);
xticklabels(xtickNames);
xtickangle(45)

% save fig       
fig = gcf;
inhFigPath = fullfile(save_path,'CM_inhibition_short_laser_catch_for_pre_inhibition.pdf');
print(inhFigPath,'-dpdf','-painters','-loose'); 

%% Summary of number of mice and sessions 
load('E:\CM_Photoinhibition_Analysis\Cross_modal_task\Performance_control\photoInh_behav')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
for i =1:length(inhSites)+1
    if i <length(inhSites)+1
        inhSite = inhSites{i};
        u = photoInh_all(strcmp(photoInh_all.inhSite, inhSite), :);
    else
        inhSite = 'all';
        u = photoInh_all;
    end
    Summary{i,1} = inhSite;
    Summary{i,2} = length(unique(u.MouseName)); % number of mice
    Summary{i,3} = size(u,1); % number of sessions
end
    
    
