function getBandICAfeature2plot(subs, conds, conds_struct, icaName, icId, groupOn, normalizeOn, checkOn, datatype, weightsOnEEG, bcorr)
% getBandICAfeature2plot(subs, conds, cond_struct, icaName, icId, groupOn, normalizeOn = 1, checkOn = 1, datatype = 'power', weightsOnEEG = 0, bcorr = '%')
% compute the average power of the specified component (by icId) for baseline and after-stimlus-onset in the specified ICA structure (icaName)
% Data before computing can be normalized (normalizedOn = 1) or not (default) in
% each trial 
% Component can be computed on raw EEG signal (weightsOnEEG = 1) or
% the extracted information (default)
% Load individual (groupOn == 0) or group IC structure (groupOn = 1) 
% Output data structure, which can be further called by plotBandICAfeature

% default
if nargin < 7 || isempty(normalizeOn) || ~ismember(normalizeOn, [0 1])
    normalizeOn = 0;
end

if nargin < 8 || isempty(checkOn) || ~ismember(checkOn, [0 1])
    checkOn = 1;
end

if nargin < 9 || isempty(datatype) || ~ ismember(datatype, {'oscillation', 'power'})
    datatype = 'power';
end

if nargin < 10 || isempty(weightsOnEEG) || ~ ismember(weightsOnEEG, [0 1])
    weightsOnEEG = 0;
end

if nargin < 11 || isempty(bcorr) || ~ ismember(bcorr, {'%', 'z', '0'}) 
    bcorr = '%';
end
bcorr_input = bcorr;
if strcmp(bcorr, '0')
    bcorr = 0;
end

% path
f_main = fileparts(which('mind_wandering2'));
f_ica = fullfile(f_main, 'ica_matfile');
cd(f_main)

% group PC
if groupOn == 1  

    % load 
    disp('Load OVALLALL ICA data ')
    load(fullfile(f_ica, ['group_', num2str(min(subs)), '_', num2str(max(subs)),'.mat']), icaName)
    eval(['ic_struct = ', icaName, ';'])
    
    weights = ic_struct.icaweights(icId,:);
    sphere = ic_struct.icasphere;
    %winv   = ic_struct.icawinv(:,icId);

    bandrange = ic_struct.bandrange;
    chans = ic_struct.chans;
end
    
for si = 1:length(subs)
    sub = subs(si);
    
    if groupOn == 0
        disp(['Load ICA data ', upper(icaName), ' of PARTICIPANT ', num2str(sub)])
        load(fullfile(f_ica, [num2str(sub), '.mat']), icaName)
        eval(['ic_struct = ', icaName, ';'])
        
        weights = ic_struct.icaweights(icId,:);
        sphere  = ic_struct.icasphere;
        %winv = ic_struct.icawinv(:, icId);
        
        bandrange = ic_struct.bandrange;
        chans = ic_struct.chans;
    end

    % compute IC per cond 
    for condi = 1:length(conds)
        cond = conds{condi};

        % load
        disp(['Load data of condition ', upper(cond)]);
        eeg = getTrials(sub, cond, conds_struct, [], chans);

        % pars
        if si == 1 && condi == 1
            srate = round(1000*length(eeg.times) / (eeg.times(end) - eeg.times(1)));       
            keyPnts = dsearchn(eeg.times', [ic_struct.times(1) -200 ic_struct.times(end)]');
        end
        
        % if no data loaded
        if size(eeg.data, 3) == 0
            disp(['Zero trials for subject ', num2str(sub), ' condition ', upper(cond)])
            continue
        end
        
        % get data matrix
        if weightsOnEEG 
            disp(['Compute ', upper(icaName), ' overall IC', num2str(icId), ' on EEG'])
            data2trans = weights * sphere * eeg.data(:,:);
            data2trans = reshape(data2trans, 1, size(eeg.data,2), size(eeg.data,3));  % nChan(=1) x nPnt x nTrial
        else 
            data2trans = eeg.data;  % nChan x nPnt x nTrial
        end
        
        % filter & hilbert-transform
        disp('Hilbert transform...')
        if si == 1 && condi == 1
            dataFilt = hilbertFilter(data2trans, srate, bandrange, checkOn);
        else 
            dataFilt = hilbertFilter(data2trans, srate, bandrange, 0);
        end

        % extract information
        if strcmpi(datatype, 'power')   
            disp('Compute raw power...')
            powerAllChans = compute_power(dataFilt, keyPnts(1):keyPnts(2), bcorr);  % baseline correction: 'db', '%', 0. 
        elseif strcmpi(datatype, 'oscillation')
            disp('Compute oscillation activity...')
            powerAllChans = real(dataFilt);  % var name keeps the same to match with prevous version when there is only "power" option
        end

        % compute IC
        if ~weightsOnEEG 
            disp(['Compute ', upper(icaName), ' overall IC', num2str(icId), ' on extracted information'])
            data_ic = weights * sphere * powerAllChans(:,:);
            data_ic = reshape(data_ic, size(powerAllChans,2), size(powerAllChans,3));
        else
            data_ic = squeeze(powerAllChans);  % nPnt x nTrial
        end

        % normalize: only for the subset part
        if normalizeOn == 1  
            data_ic(keyPnts(1):keyPnts(end),:) = normalizeK(data_ic(keyPnts(1):keyPnts(end),:), 'z');
        end

        % subset
        disp(['Subset between ', num2str(round(eeg.times(keyPnts(1)))), ' ms to ', num2str(round(eeg.times(keyPnts(end)))), ' ms'])
        data_ic = data_ic(keyPnts(1):keyPnts(end), :);  % nPnt x nTrial
        
        % save to plot
        if si == 1 && condi == 1
            data2plot = zeros(length(conds), length(subs), size(data_ic,1)); 
        end
        data2plot(condi, si, :) = mean(data_ic, 2);
        
    end  % loop over conds

end  % loop over subs

% register
dt.data = data2plot;
dt.conds = conds;
dt.condStruct = conds_struct;
dt.times = ic_struct.times;
dt.subs = subs;

% output
varName = strcat(lower(conds{1}), '_', lower(conds{end}));
f_path = fullfile(f_main, 'feat2plot_matfile');
f_temp = fullfile(f_path, [icaName, ifelse(groupOn, '_overall_', '_'), num2str(icId), ifelse(normalizeOn, '_tNorm', ''), '_bc', bcorr_input, '.mat']);
disp(['Save as data in ', f_temp])
eval([varName,' = dt;'])
if exist(f_path, 'dir') ~= 7
    disp('Create the path')
    mkdir(f_path)
end           
if exist(f_temp,'file') == 2
    disp('Overwrite the file')
    save(f_temp, varName, '-append')
else
    disp('Create the file')
    save(f_temp, varName)
end
disp('Saved')

end

