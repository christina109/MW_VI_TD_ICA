function plotIC(subs, icaName, icId, groupOn, normalizeOn, saveOn, conds2plot)
% plotIC(subs, pcaName, pcId, groupOn, normalizeOn=0, saveOn=0, conds2plot = 'all')
% plot the specified component (by icId) in the specified ICA structure (icaName)
% 1 when groupOn == 0: reading individual ICs
%   1.1 as topoplot in both average and variance
%   1.2 as waveform per individual, averaged across conditions
% 2 when groupOn == 1: reading group (overall) ICs
%   2.1 as topoplot 
%   2.2 as waveform per individual and group average, differed by
%   conditions
% Data plotted can be normalized (normalizedOn = 1) or not (default)
% Plots can be save to result path (f_result) (default not)

% default
if nargin < 5 || isempty(normalizeOn) || ~ismember(normalizeOn, [0 1])
    normalizeOn = 0;
end

if nargin < 6 || ~ismember(saveOn, [0 1]) || isempty(saveOn)
    saveOn = 0;
end

if nargin <7 || isempty(conds2plot) 
    conds2plot = {};
end

% path
f_main = fileparts(which('mind_wandering2'));
f_ica = fullfile(f_main, 'ica_matfile');
idcs   = strfind(f_main,filesep);
f_result = fullfile(f_main(1:idcs(end)-1), '4result2');
cd(f_main)

% pars
[nrow, ncol] = getPlotArrangement(length(subs));
linewidth_sub = 1.2;
shade_alpha = 0.3;
commonClimOn = 0; 

% group PC
if groupOn == 1  

    % load 
    disp('Load OVALLALL ICA data ')
    load(fullfile(f_ica, ['group_', num2str(min(subs)), '_', num2str(max(subs)),'.mat']), icaName)
    eval(['ic_struct = ', icaName, ';'])
    
    weight = ic_struct.icaweights(icId,:);
    sphere = ic_struct.icasphere;
    winv   = ic_struct.icawinv(:, icId);
    conds  = ic_struct.condcell;
    if ~isempty(conds2plot)
        [~, conds2plot] = ismember(conds2plot, conds);
    else 
        conds2plot = 1:length(conds);
    end
    
    for si = 1:length(subs)
        sub = subs(si);

        % compute & plot IC per cond 
        ylim_vals = [];
        for condi = 1:length(conds2plot)
            condid = conds2plot(condi);
            cond = conds{condid};

            % compute pc
            data = ic_struct.data{si, condid};  % nChan x nPnt x nTrial
            data_ic = weight*sphere*data(:,:);
            data_ic = reshape(data_ic, size(data,2), size(data,3));  % nPnt x nTrial

            % normalize
            if normalizeOn == 1
                data_ic = normalizeK(data_ic, 'z');
            end
            
            % plot ic wave: idv 
            % intialize
            if si == 1 && condi == 1
                fig1 = figure('position', [10 10 2000 1600]);
                set(fig1, 'PaperUnits', 'centimeters', 'PaperSize', [40, 30]);           
            end

            % create handle
            if condi == 1
                eval(['ax', num2str(si), '= subplot(nrow, ncol, ', num2str(si), ', ''parent'', fig1);'])
            end

            % prepare 
            data_mean = mean(data_ic, 2); 
            data_se = std(data_ic, [], 2)./sqrt(size(data_ic, 2)-1); 
            uE = data_mean + data_se;
            lE = data_mean - data_se;
            data_color = ic_struct.condParStruct(condid).linecolor;
            
            % save idv mean for group avg
            if si == 1 && condi == 1  % intialize
                data_ic_all = zeros([length(subs), length(conds2plot), size(data_ic,1)]);
            end
            data_ic_all(si, condi, :) = data_mean;

            % plot 
            eval(['plot(ax', num2str(si), ', ic_struct.times, data_mean, ''-'', ''linewidth'', linewidth_sub, ''color'', data_color)'])
            eval(['hold(ax', num2str(si), ', ''on'')'])
            if eval(['gca == ax', num2str(si)])
                patch([ic_struct.times, fliplr(ic_struct.times)], [lE' fliplr(uE')], data_color,'edgecolor', 'none', 'facealpha', shade_alpha)
            end

            % save for ylim
            ylim_vals = [ylim_vals, min(lE), max(uE)];

        end  % loop over conds

        % legend/axis/title
        title(num2str(sub))
        ylim(minmax(ylim_vals))
        if si == 1 
            xlabel('Time (ms)')
            ylabel('Power')
            %legend(reshape(repmat(conds(conds2plot), 2,1), length(conds2plot)*2, 1))
        end

    end  % loop over subs(1/2)

    % super graph settings
    suptitle(['overall ', upper(strrep(icaName, '_ica', '')), ' IC', num2str(icId)])

    % output
    if saveOn
        if length(conds2plot) ~= 2
            f_fig = ['overall_', strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_wave', ifelse(normalizeOn, '_norm', ''), '_idv'];
        else
            f_fig = ['overall_', strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_wave', sprintf('_%s', conds{conds2plot}), ifelse(normalizeOn, '_norm', ''), '_idv'];
        end
        print(fig1, fullfile(f_result, f_fig), '-dtiffn', '-r300')
    end

    % pc wave: group mean
    fig2 = figure('position', [10 10 2000 1600]);
    set(fig2, 'PaperUnits', 'centimeters', 'PaperSize', [40, 30]);
    
    for condi = 1:length(conds2plot)
        condid = conds2plot(condi);
        cond = conds{condid};
        
        % create handle
        ax1 = subplot(1,1,1, 'parent', fig2);
        if condi == 1
            ax1 = subplot(1,1,1, 'parent', fig2);
        end

        % prepare 
        data_mean = squeeze(mean(data_ic_all(:,condi,:), 1));
        data_se = std(data_ic_all(:,condi,:), [], 1)./sqrt(size(data_ic_all, 1)-1); data_se = squeeze(data_se);
        uE = data_mean + data_se;
        lE = data_mean - data_se;
        data_color = ic_struct.condParStruct(condid).linecolor;

        % plot 
        plot(ax1, ic_struct.times, data_mean, '-', 'linewidth', linewidth_sub, 'color', data_color)
        hold(ax1, 'on')
        patch([ic_struct.times, fliplr(ic_struct.times)], [lE' fliplr(uE')], data_color,'edgecolor', 'none', 'facealpha', shade_alpha)

        % save for ylim
        ylim_vals = [ylim_vals, min(lE), max(uE)];
        
    end  % loop over conds
    
    % legend/axis/title
    title(['overall ', upper(strrep(icaName, '_ica', '')), ' IC', num2str(icId), ' N=', num2str(length(subs))])
    ylim(minmax(ylim_vals))
    xlabel('Time (ms)')
    ylabel('Power')
    legend(reshape(repmat(conds(conds2plot), 2,1), length(conds2plot)*2, 1))
    
    % save 
    if saveOn 
        if length(conds2plot) ~= 2
            f_fig = ['overall_', strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_wave',ifelse(normalizeOn, '_norm', '')];
        else
            f_fig = ['overall_', strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_wave',sprintf('_%s', conds{conds2plot}), ifelse(normalizeOn, '_norm', '')];
        end
        print(fig2, fullfile(f_result, f_fig), '-dtiffn', '-r300')
    end

    % topoplot
    load('pars_EEG.mat', 'chanlocs')

    fig3 = figure('position', [10 10 2000 1600]);
    set(fig3, 'PaperUnits', 'centimeters', 'PaperSize', [40, 30]);
    topoplot(winv, chanlocs, 'electrodes', 'labels');
    colorbar
    title(['overall ', upper(strrep(icaName, '_ica', '')), ' IC', num2str(icId),' N=', num2str(length(subs))])

    % save 
    if saveOn 
        f_fig = ['overall_', strrep(icaName, '_ica', ''), '_ic', num2str(icId)];
        print(fig3, fullfile(f_result, f_fig), '-dtiffn', '-r300')
    end

%%
% individual IC
else  
    
    for si = 1:length(subs)
        sub = subs(si);

        % load 
        disp(['Load ICA data ', upper(icaName), ' of PARTICIPANT ', num2str(sub)])
        load(fullfile(f_ica, [num2str(sub), '.mat']), icaName)
        eval(['ic_struct = ', icaName, ';'])

        % intialize
        if si == 1
            winvs = zeros(length(subs), size(ic_struct.icawinv,1)); 
            conds = ic_struct.condcell;
            if ~isempty(conds2plot)
                [~, conds2plot] = ismember(conds2plot, conds);
            else
                conds2plot = 1:length(conds);
            end
        end
        
        % register
        winvs(si, :) = ic_struct.icawinv(:, icId);
        loadings = ic_struct.icaweights(icId,:)*ic_struct.icasphere;

        % compute & plot PC per cond 
        ylim_vals = [];
        for condi = 1:length(conds2plot)
            condid = conds2plot(condi);
            cond = conds{condid};

            % compute pc
            data = ic_struct.data{condid};  % nChan x nPnt x nTrial
            data_ic = loadings*data(:,:);
            data_ic = reshape(data_ic, size(data,2), size(data,3));  % nPnt x nTrial

            % normalize
            if normalizeOn == 1
                data_ic = normalizeK(data_ic, 'z');
            end

            % plot pc wave
            % intialize
            if si == 1 && condi == 1
                fig1 = figure('position', [10 10 2000 1600]);
                set(fig1, 'PaperUnits', 'centimeters', 'PaperSize', [40, 30]);
            end

            % create handle
            if condi == 1
                eval(['ax', num2str(si), '= subplot(nrow, ncol, ', num2str(si), ', ''parent'', fig1);'])
            end

            % prepare 
            data_mean = mean(data_ic, 2); 
            data_se = std(data_ic, [], 2)./sqrt(size(data_ic, 2)-1); 
            uE = data_mean + data_se;
            lE = data_mean - data_se;
            data_color = ic_struct.condParStruct(condid).linecolor;

            % plot 
            eval(['plot(ax', num2str(si), ', ic_struct.times, data_mean, ''-'', ''linewidth'', linewidth_sub, ''color'', data_color)'])
            eval(['hold(ax', num2str(si), ', ''on'')'])
            if eval(['gca == ax', num2str(si)])
                patch([ic_struct.times, fliplr(ic_struct.times)], [lE' fliplr(uE')], data_color,'edgecolor', 'none', 'facealpha', shade_alpha)
            end

            % save for ylim
            ylim_vals = [ylim_vals, min(lE), max(uE)];

        end  % loop over conds

        % legend/axis/title
        title(num2str(sub))
        ylim(minmax(ylim_vals))
        if si == 1 
            xlabel('Time (ms)')
            ylabel('Power')
            %legend(reshape(repmat(conds, 2,1), length(conds)*2, 1))
        end

    end  % loop over subs(1/2)

    % super graph settings
    suptitle([upper(strrep(icaName, '_ica', '')), ' IC', num2str(icId)])

    % output
    if saveOn
        if length(conds2plot) ~= 2
            f_fig = [strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_wave', ifelse(normalizeOn, '_norm', ''), '_idv'];
        else
            f_fig = [strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_wave', sprintf('_%s', conds{conds2plot}), ifelse(normalizeOn, '_norm', ''), '_idv'];
        end
        print(fig1, fullfile(f_result, f_fig), '-dtiffn', '-r300')
    end

    % topoplot: individual
    load('pars_EEG.mat', 'chanlocs')
    fig2 = figure('position', [10 10 2000 1600]);
    set(fig2, 'PaperUnits', 'centimeters', 'PaperSize', [40, 30]);
    clim_rng = minmax(reshape(winvs, 1, numel(winvs)));

    for si = 1:length(subs)
        subplot(nrow, ncol, si, 'parent', fig2)

        if commonClimOn 
            topoplot(winvs(si,:), chanlocs, 'maplimits', clim_rng);
        else
            topoplot(winvs(si,:), chanlocs);
        end

        % graph settings
        title([num2str(subs(si))])
        if commonClimOn 
            if si == 1
                colorbar
            end
        %else
        %    colorbar
        end

    end  % loop over subs(2/2)

    % super graph settings
    suptitle([upper(strrep(icaName, '_ica', '')), ' IC', num2str(icId)])

    % save 
    if saveOn 
        f_fig = [strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_idv'];
        print(fig2, fullfile(f_result, f_fig), '-dtiffn', '-r300')
    end

    % topoplot: group
    % norm weight per individual
    winvs_norm = normalizeK(winvs, 'range', 1);

    fig3 = figure('position', [10 10 2000 1600]);
    set(fig3, 'PaperUnits', 'centimeters', 'PaperSize', [40, 30]);

    subplot(1,2,1)
    topoplot(mean(winvs_norm,1), chanlocs);
    colorbar
    title(['normalized ', upper(strrep(icaName, '_ica', '')), ' IC', num2str(icId)])

    subplot(1,2,2)
    topoplot(std(winvs_norm,[],1), chanlocs);
    colorbar
    title(['std N=', num2str(length(subs))])

    % save 
    if saveOn 
        f_fig = [strrep(icaName, '_ica', ''), '_ic', num2str(icId), '_norm'];
        print(fig3, fullfile(f_result, f_fig), '-dtiffn', '-r300')
    end
end

end