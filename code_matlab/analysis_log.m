
% main working directory
f_main = fileparts(which('mind_wandering2'));
f_result = fileparts(which('mind_wandering2_results'));
cd(f_main)

%% Jul 19, 2019
% prepare and plot IC feature behavior

% ***   prepare the data   ***
subs = 301:330;
trialNormOn = 0;

conds = {'eoa', 'ioa'};
load('pars_conds', 'conds_raw')
%ics =  setdiff([ics_t, ics_tv, ics_tm, ics_tvm], ics2ex);  % vars from Jun 20
ics =  setdiff([ics_t, ics_tm], ics2ex); 
for ici = 1:length(ics)
    %getBandICAfeature2plot(subs, conds, conds_raw, 'alpha_ica_oscillation', ics(ici), 1, trialNormOn, 0, 'power', 1)
    getBandICAfeature2plot(subs, conds, conds_raw, 'alpha_ica_oscillation', ics(ici), 1, trialNormOn, 0, 'power', 1, '0')
end

conds = {'nf', 'f'};
%conds = {'ot', 'mw'};
load('pars_conds', 'conds_beh')
%ics = setdiff([ics_v, ics_tv, ics_tvm], ics2ex);
%ics = setdiff([ics_m, ics_tm, ics_tvm], ics2ex);
ics = setdiff([ics_v], ics2ex);
%ics = setdiff([ics_m, ics_tm], ics2ex);
for ici = 1:length(ics)
    %getBandICAfeature2plot(subs, conds, conds_beh, 'alpha_ica_oscillation', ics(ici), 1, trialNormOn, 0, 'power', 1)
    getBandICAfeature2plot(subs, conds, conds_beh, 'alpha_ica_oscillation', ics(ici), 1, trialNormOn, 0, 'power', 1,'0')
end

% ***   plot the data   ***
linecolors = {{[30,144,255]./255, 'r'}, {'b', [255,20,147]./255}, {[50,205,50]./255, [255,140,0]./255}};

condsets = {{'eoa', 'ioa'}, {'nf', 'f'}, {'ot', 'mw'}};

%icsets = {ics_t, ics_v, ics_m, ics_tv, ics_tm, ics_tvm};
%icconds = {1, 2, 3, [1 2], [1 3], 1:3};
icsets = {ics_t, ics_v, ics_m, ics_tm};
icconds = {1, 2, 3, [1 3]};
nrow = 4; %nrow = 4;
ncol = 5; %ncol = 6;

figure;
% spi = 0;
for seti = 1:length(icsets)
    icset = icsets{seti};
    iccond = icconds{seti};
    ics =  setdiff(icset, ics2ex); 
    spi = 0;

    for ici = 1:length(ics)
        ic = ics(ici);
        %f_load = fullfile('feat2plot_matfile', ['alpha_ica_oscillation_overall_', num2str(ic), '.mat']);
        f_load = fullfile('feat2plot_matfile', ['alpha_ica_oscillation_overall_', num2str(ic), '_bc0.mat']);

        for condsi = 1:length(iccond)
            conds = condsets{iccond(condsi)};
            colors = linecolors{iccond(condsi)};
            load(f_load, [conds{1}, '_', conds{2}])
            eval(['data = ', [conds{1}, '_', conds{2}], ';'])

            % scan for "all zeros" 
            % some participants lack one of the conditions due to
            % self-reports
            temp = mean(data.data, 3); 
            subsin = and(temp(1,:), temp(2,:));  % logical vector containing "non-zero" subjects

            % prepare 
            y = mean(data.data(:, subsin, :), 2); y = squeeze(y);
            se = std(data.data(:, subsin, :),[], 2)./sqrt(sum(subsin)-1); se = squeeze(se);
            ymin = y - se;
            ymax = y + se;
            yrng = [min(min(ymin)), max(max(ymax))];

            % plot
            spi = spi + 1; 
            subplot(nrow, ncol, ncol*(seti-1)+spi)
            for condi = 1:2
                patch([data.times, flip(data.times)], [ymin(condi,:), flip(ymax(condi,:))], colors{condi}, 'edgecolor', 'none', 'facealpha', 0.3)
                hold on
                plot(data.times, y(condi,:), 'color', colors{condi}, 'linewidth', 1)
                title(['IC' num2str(ic)])
            end

            % axis settings
            %ylim(yrng)
            xlim([-600 600])
            %set(gca, 'yTick', [0 round(yrng(2)/2) floor(yrng(2))])
            %set(gca, 'yTick', [ceil(yrng(1)) round((yrng(2)-yrng(1))/2+yrng(1)) floor(yrng(2))])
            ylim([floor(yrng(1)) ceil(yrng(2))])
            set(gca, 'yTick', [floor(yrng(1)) ceil(yrng(2))])
            
            % stimlus onset line
            plot([0 0], [floor(yrng(1)) ceil(yrng(2))], 'color', [0.3 0.3 0.3])
            %plot([-600 600], [0 0], 'color', [0.3 0.3 0.3])
            
            %if spi == 1
            if seti == 1 && spi == 1
                %ylabel('Change in percentage')
                ylabel('Raw power') % check the unit
            end
            %if spi <= (nrow-1)*ncol
            if seti < nrow
                set(gca,'XTick',[])
            end
            %if spi == (nrow-1)*ncol +1
            if seti == nrow && spi == 1
                xlabel('Time (ms)')
            end
        end  % loop over condi
    end  % loop over ici
end  % loop over seti

% deal with missing values
%% Jun 20, 2019
% run dipole fitting with overall ICs
eeglab

icaName = 'alpha_ica_oscillation';
disp('Load ICA result structure...')

load(fullfile('ica_matfile', 'group_301_330.mat'), icaName);
eval(['icaStruct = ', icaName, ';'])

% path pars
eeglabPath = fileparts(which('eeglab'));
bemPath = fullfile(eeglabPath, 'plugins', 'dipfit2.3', 'standard_BEM'); 

% arrange to EEG format
disp('Initialize EEG...')
EEG = struct();
disp('Concatenate data...')
for si = 1:length(icaStruct.subs)
    for ci = 1:length(icaStruct.condcell)
        data2copy = icaStruct.data{si,ci};
        if si * ci == 1
            EEG.data = data2copy;
        else
            EEG.data = cat(3, EEG.data, data2copy);
        end
    end
end
% check
nTrial = 0;
for si = 1:length(icaStruct.subs)
    for ci = 1:length(icaStruct.condcell)
        nTrial = size(icaStruct.data{si,ci},3) + nTrial;
    end
end
sprintf('The overall data size are: %i %i %i', size(EEG.data))
if size(EEG.data,3) ~= nTrial
    error('Wrongly loaded')
end

% add other pars
disp('Add other parameters...')
EEG.icawinv = icaStruct.icawinv;
EEG.icaweights = icaStruct.icaweights;
EEG.icasphere = icaStruct.icasphere;
EEG.icaact = [];
EEG.times = icaStruct.times;
EEG.xmax = max(EEG.times)/1000;
EEG.xmin = min(EEG.times)/1000;
load('pars_EEG', 'chanlocs', 'srate', 'chaninfo')
EEG.chanlocs = chanlocs;
EEG.icachansind = 1:length(chanlocs);
EEG.nbchan = length(chanlocs);
EEG.dipfit = [];
EEG.setname = 'concatenated data for ICA';
EEG.trials = size(EEG.data,3);
EEG.pnts = size(EEG.data,2);
EEG.srate = srate; 
EEG.chaninfo = chaninfo;
disp('EEG structure is ready.')
EEG_bk = EEG;

% ***   dipole fitting   ***
% Setting up DIPFIT model and preferences
eeglab
EEG = EEG_bk;  % start eeglab will intialize EEG
disp('Set DIPFIT model and preferences...')
EEG = pop_dipfit_settings(EEG, 'hdmfile',fullfile(bemPath, 'standard_vol.mat'),'coordformat','MNI',...
                          'mrifile',fullfile(bemPath, 'standard_mri.mat'),...
                          'chanfile',fullfile(bemPath, 'elec', 'standard_1005.elc'),...
                          'chansel',1:length(EEG.chanlocs));
EEG = pop_chanedit(EEG, 'lookup',fullfile(bemPath, 'elec', 'standard_1005.elc'));

% Initial fitting - Scanning on a coarse-grained grid
%EEG = pop_dipfit_gridsearch(EEG, 1:32, -85:17:85, -85:17:85, 0:17:85, 0.4);
%pop_dipplot(EEG, 1:length(EEG.chanlocs), 'mri',fullfile(bemPath, 'standard_mri.mat'),'normlen','on');

% coarse+fine fitting 
ics2fit = 1:32;
threshold = 100;  % in % to filter out fitting results with higher residul vars
EEG = pop_multifit(EEG, ics2fit, 'threshold', threshold); 

% fine tune for IC1 
EEG = pop_dipfit_nonlinear(EEG);  % adjust to 2 dips for IC 1

% save 
save(['EEG4dipfit_', icaName,'.mat'], 'EEG')

% ***   plot dipole fitting results   ***
% load
f_dipfit = ['EEG4dipfit_', icaName,'.mat'];
load(f_dipfit, 'EEG')

% plot results
%t-tests to ideal chance level
%icsall = [1  2  3  4  5  6  7  8  9 10 11 12 15 16 17 18 19 20 21 23 24 26 28 30];
%ics2ex = [5 10 20 21 23 24 26 28 30];
%ics_tvm = 2;
%ics_tv = [3 7 15];
%ics_tm = [4 5 9];
%ics_t = [6 8 10 11 16 19 18 23 28 30];
%ics_v = [1 12 26];
%ics_m = [17 20 21 24];

icsall = [1  2  3  4  5  6  7  9 15 17 26];
ics2ex = [5 10 20 21 23 24 26 28 30];
ics_tm = [2 4];
ics_t = [3 5 6 7 9 15];
ics_v = [1 26];
ics_m = [17];

ics2plot = setdiff(ics_tm, ics2ex);
%ics2plot = icsall;
rvrange = [0 40];  
%pop_dipplot(EEG, ics2plot ,'rvrange', rvrange, 'mri', fullfile(bemPath, 'standard_mri.mat'),...
%            'cornermri','on','axistight','on','projimg','on','projlines','on','pointout','on','normlen','on','view',[51 18] );
pop_dipplot(EEG, ics2plot ,'rvrange', rvrange, 'mri', fullfile(bemPath, 'standard_mri.mat'),...
            'summary','on','num','on','normlen','on');  % filter by rvrange not work?

t_fig = 'Dipole fitting with ICs';
[nrow,ncol] =  getPlotArrangement(length(ics2plot));
pop_topoplot(EEG,0, ics2plot ,t_fig, [nrow,ncol], 1, 'electrodes','on');
% pop_topoplot(EEG,0, ics2plot ,t_fig, [4 6], 1, 'electrodes','on');

% note: do fine fitting if there is a need 

% get coordinates
mat = zeros(length(icsall)+1, 7);
for ici = 1:length(icsall)
    ic = icsall(ici);
    if ic == 1
        mat(ici:ici+1,1) = [ic ic];
        mat(ici:ici+1,2:4) = round(EEG.dipfit.model(ic).posxyz);
        %mat(ici:ici+1,5:7) = round(EEG.dipfit.model(ic).momxyz);
    else 
        mat(ici+1, 1) = ic;
        mat(ici+1,2:4) = round(EEG.dipfit.model(ic).posxyz(1,:));
        %mat(ici+1,5:7) = round(EEG.dipfit.model(ic).momxyz(1,:));
    end
end

% output 
csvwrite(fullfile(f_result, 'dip_coords.csv'), mat)

% convert to tal
temp = mat(:,2:4);
temp = mni2tal(temp);
mat_tal = [mat(:,1) round(temp)];
% cuixuFindStructure()

