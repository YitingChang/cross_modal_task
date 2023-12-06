%%% behavioral performance
% Input: MSessionExplorer of all inhibition sessions for each mouse
% Output: Session-by-session summary of behavioral performance (lick info) 
% For each session, behavioral performance must pass the following criteria.
            % Criteria: 
            % Hit rate: lick probability of hit trials
            % Catch rate: lick probability of catch trials (no stimulus and laser)
            % Laser catch rate: lick probability of laser catch trials
            
            % Without inhibition, 
            % 1. hit rate >= 35% 
            % 2. overall performance >= 65%
            % Without stimuli,
            % 3. laser catch rate (short and long) <= hit rate-5%
            

% Directory
Mice = {'YT091' 'YT092' 'YT093' 'YT094' 'YT095'};
% Disks = {'D:\' 'E:\' 'E:\' 'F:\' 'F:\' 'E:\' 'E:\'};

for m=1:length(Mice)
    MouseName = Mice{m};
%     MainDir = Disks{m};
    seDir = ['F:\', MouseName, '\MSessionExplorer_inhibition'];
%     seFilePaths = MBrowse.Files(seDir);
    cd(seDir)
    seFileLists = struct2cell(dir('*se*'));
    seFileNames = seFileLists(1,:);
    seFilePaths = cellfun(@(x) fullfile(seDir, x),seFileNames,'UniformOutput',false);
    behavSummary = CM.Inhibition.GetBehavInfo_tDetection(seFilePaths);

    % Save behavior summary
    photoInhPaths = fullfile('E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control', ...
        [MouseName,'_photoinhibition']);
    save(photoInhPaths, 'behavSummary');
end

%% Concatenate behavior summary for all mice 
photoInhDir = ['E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control'];
photoInhPaths = MBrowse.Files(photoInhDir);

photoInh_all = table();
for i = 1: length(photoInhPaths)
    load(photoInhPaths{i})
    % Sessions with qualified perfromance are selected
    isPassed = behavSummary.isPassed;
    photoInh_all = [photoInh_all; behavSummary(isPassed,:)];
end
save('E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control\photoInh_behav', 'photoInh_all');

%% Lick probability for individual mice
% Select lick responses based on brain regions of interest
load('E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control\photoInh_behav')
Genotypes = {'PV-Cre-Ai32'}; 
% Genotypes = {'PV-Cre-Ai32' 'VGAT-ChR2' 'all'};
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control\Figures_lick_prob\';
S=5; % tactile, no stimulus
O=3; % no, pre-, post- inhibition 
nboot = 10000; % number for bootstrapping
titles_tt = {'Tactile trials' 'Catch trials'};
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
            % Bootstrap lick responses for all trial types
            for s=1:S
                for o=1:O
                    v=[];
                    % Concatenation
                    for i = 1:height(m)
                        new_v = m{i,5}{1,1}{s,o}; 
                        v = [v; new_v];
                    end
                    if length(v) >1
                        [ci bootstat]= bootci(nboot,@mean,v);  
                        MEAN = mean(bootstat);
                        data_bootstrapped(s,o,1) = MEAN;
                        data_bootstrapped(s,o,2) = ci(1);
                        data_bootstrapped(s,o,3) = ci(2);
                    else
                        data_bootstrapped(s,o,1) = NaN;
                        data_bootstrapped(s,o,2) = NaN;
                        data_bootstrapped(s,o,3) = NaN;
                    end
                end  
             end
            
            % plot for individual mouse
            figure('Position', [0,0, 600, 300]); % set figure size
            %%%%% tactile and catch trials %%%%%
            for s=1:2
                subplot(1,3,s+1)
                for o=1:O
                    x = o;
                    y1 = data_bootstrapped(s,o,1);
                    y2 = data_bootstrapped(s,o,2);
                    y3 = data_bootstrapped(s,o,3);
                    errorbar(x,y1,(y1-y2),(y3-y1),'o','Color','b',...
                    'MarkerSize',6,'MarkerFaceColor','b', 'MarkerEdgeColor', 'b'); hold on
                end
                set(gca, 'box','off','TickDir','out');
                ylim([0 1]);
                ax = gca; ax.YAxis.Visible = 'off';
                xlim([0.5 3.5]);  
                xticklabels({'no manipulation' 'pre-inhibition' 'post-inhibition'});
                xtickangle(45)
                title(titles_tt{s})
            end
            
            %%%%% laser catch trial %%%%%
            subplot(1,3,1)
            for s=3:5
                for o= [1,3]
                    if o==1
                        x = s-3;
                    else
                        x = s;
                    end
                    y1 = data_bootstrapped(s,o,1);
                    y2 = data_bootstrapped(s,o,2);
                    y3 = data_bootstrapped(s,o,3);
                    errorbar(x,y1,(y1-y2),(y3-y1),'o','Color','b',...
                    'MarkerSize',6,'MarkerFaceColor','b', 'MarkerEdgeColor', 'b')
                    hold on
                end
            end
            set(gca, 'box','off','TickDir','out');
            ylim([0 1]);
            ylabel('lick probability');
            xlim([-0.5 5.5]); 
            xticks(0:5)
            xticklabels({'all\_ITI' 'short\_ITI' 'long\_ITI' 'all\_laser only' 'short\_laser only' 'long\_laser only'})
            xtickangle(45)
            title('Laser catch trials')
              
            sgtitle([mouseName '\_' inhSite], 'FontSize',14);
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

load('E:\CM_Photoinhibition_Analysis\Tactile_detection\photoInh_behav')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 

Genotypes = {'PV-Cre-Ai32'}; 
S=5; % tactile, no stimulus
O=3; % no, pre-, post- inhibition 
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
                    % Tactile and catch trials 
                    for s=1:2
                        for o=1:O
                            % resample trials
                            lick_tt = v{s,o};
                            resampled_lick_tt = randsample(lick_tt, length(lick_tt), true); 
                            lick_prob(s,o) = mean(resampled_lick_tt); % lick probability of each trial type
                            if o > 1
                                delta_lick_prob_session(session,s,o-1) = lick_prob(s,o) - lick_prob(s,1);
                                % difference in lick probability (inhibition - no inhibition)
                            end
                        end
                    end
                    % Laser catch trials
                    for s=3:5
                        for o= [1,3]
                            % resample trials
                            lick_tt = v{s,o};
                            resampled_lick_tt = randsample(lick_tt, length(lick_tt), true); 
                            lick_prob(s,o) = mean(resampled_lick_tt); % lick probability of each trial type
                            if o > 1
                                delta_lick_prob_session(session,s,1) = lick_prob(s,o) - lick_prob(s,1);
                                % difference in lick probability (inhibition - no inhibition)
                            end
                        end
                    end
                    
                    % Laser short catch trials (iti vs. response window after laser offset)
                    delta_lick_prob_session(session,4,2) = lick_prob(2,2) - lick_prob(4,1);
                    clear lick_prob
                end
                delta_lick_prob_mouse(mouse,1:S,1:(O-1)) = mean(delta_lick_prob_session,1); % average across sessions
                clear delta_lick_prob_session
            end
            delta_lick_prob(n,1:S,1:(O-1)) = mean(delta_lick_prob_mouse,1); % average across mice
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
photoInhPaths = fullfile('E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control', ... 
    'photoinhibition_bootstrapped_PV-Cre-Ai32');

save(photoInhPaths, 'photoInh_bootstrapped'); 
% clear photoInh_bootstrapped 

%% Plotting
load('E:\CM_Photoinhibition_Analysis\Tactile_detection\photoinhibition_bootstrapped_PV-Cre-Ai32')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
Genotypes = {'PV-Cre-Ai32'};
titles_tt = {'pre-inhibition' 'post-inhibition'};
ylabels_tt = {{'Tactile trials' 'delta right lick prob'} {'Catch trials' 'delta right lick prob'}};
ylabels_laser = {{'All laser trials' 'delta right lick prob'}...
    {'Short laser trials' 'delta right lick prob'} {'Long laser trials' 'delta right lick prob'}};
colors_tt = {'b' 'b'};
xtickNames = {'sham' 'S1/S2' 'MM' 'ALM' 'AMM' 'Prt'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Tactile_detection\Figures_delta_lick_prob';
S=5; % tactile, no stimulus
O=3; % no, pre-, post- inhibition 
% catch trial: all (s=3), short (s=4), long (s=5)
% o=1: ITI of laser catch trials; o=3: respond window of laser catch trials
nboot = 10000; % number for bootstrapping

% Tactile and catch trials
figure('Position', [0,0, 600, 600]); % set figure size
for s=1:2   
    for o=1:O-1
        subplot(2,2,(s-1)*(O-1)+o)
        for e = 1:length(inhSites) 
            isInhSites = cellfun(@(x) strcmp(x,inhSites{e}),photoInh_bootstrapped(:,2));
            delta_lick_prob_bootstrapped_tt = photoInh_bootstrapped{isInhSites,3}(:,s,o);
            sorted_delta_lick_prob = sort(delta_lick_prob_bootstrapped_tt);   
            y1 = mean(sorted_delta_lick_prob); 
            y2 = sorted_delta_lick_prob(nboot*0.025); 
            y3 = sorted_delta_lick_prob(nboot*0.975); 
            x = e;
            errorbar(x,y1,(y1-y2),(y3-y1),'o','Color',colors_tt{s},...
            'MarkerSize',6,'MarkerFaceColor',colors_tt{s}, 'MarkerEdgeColor',colors_tt{s})
            hold on
            CI(o,e,1) = y2;
            CI(o,e,2) = y3;
        end
        set(gca, 'box','off','TickDir','out');
        ylim([-1 0.5]);
        set(gca,'Ytick',-1:0.5:0.5); 
        ylabel(ylabels_tt{s});
        yline(0, '--k');
        xlim([0.5 (length(inhSites)+0.5)]);  
        xticks([1:length(inhSites)]);
        xticklabels(xtickNames);
        xtickangle(45)
        title(titles_tt{o})
    end
end
% save fig       
fig = gcf;
inhFigPath = fullfile(save_path,'tacDetection_inhibition_tactile_and_catch.pdf');
print(inhFigPath,'-dpdf','-painters','-loose');   

% Laser catch trials
figure('Position', [0,0, 300, 900]); % set figure size
for s=3:5
    subplot(3,1,s-2)
    for e = 1:length(inhSites) 
        isInhSites = cellfun(@(x) strcmp(x,inhSites{e}),photoInh_bootstrapped(:,2));
        delta_lick_prob_bootstrapped_tt = photoInh_bootstrapped{isInhSites,3}(:,s,1);
        sorted_delta_lick_prob = sort(delta_lick_prob_bootstrapped_tt);   
        y1 = mean(sorted_delta_lick_prob); 
        y2 = sorted_delta_lick_prob(nboot*0.025); 
        y3 = sorted_delta_lick_prob(nboot*0.975); 
        x = e;
        errorbar(x,y1,(y1-y2),(y3-y1),'o','Color','b',...
        'MarkerSize',6,'MarkerFaceColor','b', 'MarkerEdgeColor','b')
        hold on
    end
    set(gca, 'box','off','TickDir','out');
    ylim([-1 0.5]);
    set(gca,'Ytick',-1:0.5:0.5); 
    ylabel(ylabels_laser{s-2});
    yline(0, '--k');
    xlim([0.5 (length(inhSites)+0.5)]);  
    xticks([1:length(inhSites)]);
    xticklabels(xtickNames);
    xtickangle(45)
end
% save fig       
fig = gcf;
inhFigPath = fullfile(save_path,'tacDetection_inhibition_laser_catch.pdf');
print(inhFigPath,'-dpdf','-painters','-loose'); 

% Laser short catch trials (iti vs. response window after laser offset) Figure 5F
figure('Position', [0,0, 300, 300]); % set figure size
for e = 1:length(inhSites) 
    isInhSites = cellfun(@(x) strcmp(x,inhSites{e}),photoInh_bootstrapped(:,2));
    delta_lick_prob_bootstrapped_tt = photoInh_bootstrapped{isInhSites,3}(:,4,2);
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
ylabel('delta right lick prob');
yline(0, '--k');
xlim([0.5 (length(inhSites)+0.5)]);  
xticks([1:length(inhSites)]);
xticklabels(xtickNames);
xtickangle(45)
title({'Short laser catch trials', '(iti vs. window after laser offset)'})
% save fig       
fig = gcf;
inhFigPath = fullfile(save_path,'tacDetection_inhibition_short_laser_catch.pdf');
print(inhFigPath,'-dpdf','-painters','-loose'); 


%% sensitivity index
% Hit - catch trials

% Select lick responses based on brain regions of interest
load('E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control\photoInh_behav')
Genotypes = {'PV-Cre-Ai32'}; 
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control\Figures_sensitivity\';

S=5; % tactile, no stimulus
O=3; % no, pre-, post- inhibition 
% catch trial: all (s=3), short (s=4), long (s=5)
% o=1: ITI of catch trials; o=3: laser catch trials
optoNames = {'control' 'pre-inh.' 'post-inh'};
sgtitles = {'sham' 'S1\_S2' 'MM' 'ALM' 'AMM' 'Prt'};
ylabels = {'Hit rate' 'FA rate' 'Sensitivity'};

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
                for o=1:O
                    Hit = v{1,o}; 
                    pHit = mean(Hit);% hit rate
                    FA = v{2,o}; 
                    pFA = mean(FA);% false alarm rate
                    sensitivity = pHit - pFA;
%                         dpri = dprime(pHit,pFA,length(Hit),length(FA));
                    data_session(1,session,o) = pHit;
                    data_session(2,session,o) = pFA;
                    data_session(3,session,o) = sensitivity;
                end
            end
            data_mouse(1:3,mouse,1:O) = mean(data_session,2);
            clear data_session
        end
        % plot all mice
        figure('Position', [0,0, 300, 900]); % set figure size
        for k=1:3
            subplot(3,1,k)
            for o=1:O
                y = data_mouse(k,:,o);
                x = repmat(o,length(y),1);
                scatter(x,y,'b');hold on;
            end
            line = squeeze(data_mouse(k,:,:))';
            plot(line,'Color',[0.7 0.7 0.7])  
            xlim([0.5 3.5])
            set(gca,'Xtick',1:3);
            xticklabels(optoNames)
            ylabel(ylabels{k})
            if k == 3
                ylim([-1 1])
            else
                ylim([0 1])
            end
        end
        sgtitle(sgtitles{e})
        % Save fig
        fig = gcf;
        inhFigPath = fullfile(save_path,['sensitivity_' inhSites{e} '.pdf']);
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

clearvars -except photoInh_all
Genotypes = {'PV-Cre-Ai32'}; 
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 

S=5; % tactile, no stimulus
O=3; % no, pre-, post- inhibition 
nboot = 10000; % number for bootstrapping     
seedNum = [1 2 3 4 5 6 7];
% for sn = 4:length(seedNum)
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
            rng(8)
%             rng(seedNum(sn));
            for n = 1:nboot
                resampled_mice = randsample(mouseNames, length(mouseNames), true); 
                for mouse = 1:length(resampled_mice)
                    m = u(strcmp(u.MouseName, resampled_mice{mouse}), :);
                    % resample session
                    for session = 1:height(m)
                        v = m{randi(height(m)),5}{1,1};   
                        for o=1:O
                            % resample hit trials      
                            Hit = v{1,o}; 
                            resampled_Hit = randsample(Hit, length(Hit), true); 
                            pHit = mean(resampled_Hit);% hit rate
                            % resample FA trials
                            FA = v{2,o}; 
                            resampled_FA = randsample(FA, length(FA), true);
                            pFA = mean(FA);% false alarm rate
                            sensitivity = pHit - pFA;
    %                             dpri = dprime(pHit,pFA,length(Hit),length(FA));
                            st_session(session,o) = sensitivity; 
                            if o > 1
                                delta_st_session(session,o-1) = sensitivity - st_session(session,1);
                            end
                        end
                    end
                    delta_st_mouse(mouse,1:O-1) = mean(delta_st_session,1); 
                    clear st_session delta_st_session
                end
                delta_sensitivity(n,1:O-1) = mean(delta_st_mouse,1);  
                clear delta_st_mouse
            end
            delta_sensitivity_bootstrapped{(g-1)*length(inhSites)+e,1} = Genotype;
            delta_sensitivity_bootstrapped{(g-1)*length(inhSites)+e,2} = inhSite;
            delta_sensitivity_bootstrapped{(g-1)*length(inhSites)+e,3} = delta_sensitivity;
            clear delta_sensitivity
        end
    end

    % Save 
    photoInhPaths = fullfile('E:\CM_Photoinhibition_Analysis\Tactile_detection\Performance_control', ... 
        'delta_sensitivity_bootstrapped_PV-Cre-Ai32');
    save(photoInhPaths, 'delta_sensitivity_bootstrapped'); 
%     clear delta_sensitivity_bootstrapped
% end       
%% Figure 5C
% delta sensitivity
% analysis setting
load('E:\CM_Photoinhibition_Analysis\Tactile_detection\delta_sensitivity_bootstrapped_PV-Cre-Ai32')
inhSites = {'sham' 'S1' 'MM' 'ALM' 'AMM' 'Prt'}; 
save_path = 'E:\CM_Photoinhibition_Analysis\Tactile_detection\Figures_delta_sensitivity\';
S=5; % tactile, no stimulus
O=3; % no, pre-, post- inhibition 
nboot = 10000; % number for bootstrapping  

% figure setting
titles = {'pre-inhibition' 'post-inhibition'};
xtickNames = {'sham' 'S1/S2' 'MM' 'ALM' 'AMM' 'Prt'};
% colors_inhibition = {[0.65 0.9 0.5] [0.25 0.5 0.05]};

figure('Position', [0,0, 700, 300]); % set figure size
for o=1:O-1
    subplot(1,2,o)
    for e = 1:length(inhSites) 
        isInhSites = cellfun(@(x) strcmp(x,inhSites{e}),delta_sensitivity_bootstrapped(:,2));
        delta_sensitivity_bootstrapped_tt = delta_sensitivity_bootstrapped{isInhSites,3}(:,o);
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
%         ylim([-3.5 1]);
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
inhFigPath = fullfile(save_path,'tacDetection_inhibition_delta_sensitivity_v2.pdf');
print(inhFigPath,'-dpdf','-painters','-loose');  

%% Summary of number of mice and sessions 
load('E:\CM_Photoinhibition_Analysis\Tactile_detection\photoInh_behav')
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