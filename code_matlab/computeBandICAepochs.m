function computeBandICAepochs(subs, icaName, icId, groupOn, normalizeOn, checkOn, datatype, weightsOnEEG)
% computeBandICAepochs(subs, icaName, icId, groupOn, normalizeOn=0, checkOn = 1, datetype = 'power', weightsOnEEG = 0)
% compute the average power of the specified component (by icId) for baseline and after-stimlus-onset in the specified ICA structure (icaName)
% save as feature files 
% Data before computing can be normalized (normalizedOn = 1) or not (default) in
% each trial 
% Component can be computed on raw EEG signal (weightsOnEEG = 1) or
% the extracted information (default)
% Load individual (groupOn == 0) or group IC structure (groupOn = 1) 

% default
if nargin < 5 || isempty(normalizeOn) || ~ismember(normalizeOn, [0 1])
    normalizeOn = 0;
end

if nargin < 6 || isempty(checkOn) || ~ismember(checkOn, [0 1])
    checkOn = 1;
end

if nargin < 7 || isempty(datatype) || ~ ismember(datatype, {'oscillation', 'power'})
    datatype = 'power';
end

if nargin < 8 || isempty(weightsOnEEG) || ~ ismember(weightsOnEEG, [0 1])
    weightsOnEEG = 0;
end

% path
f_main = fileparts(which('mind_wandering2'));
f_ica = fullfile(f_main, 'ica_matfile');
f_output = fullfile(f_main, 'measure_matfile');
cd(f_main)

% pars
stimsepOn = 0;
if stimsepOn 
    wincut = [-200 0];  % [baseline end, stim onset]
else
    winwid = 100;
end

% group PC
if groupOn == 1  

    % load 
    disp('Load OVALLALL ICA data ')
    load(fullfile(f_ica, ['group_', num2str(min(subs)), '_', num2str(max(subs)),'.mat']), icaName)
    eval(['ic_struct = ', icaName, ';'])
    
    weights = ic_struct.icaweights(icId,:);
    sphere = ic_struct.icasphere;
    %winv   = ic_struct.icawinv(:,icId);
    
    conds  = ic_struct.condcell;
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
        conds = ic_struct.condcell;

        bandrange = ic_struct.bandrange;
        chans = ic_struct.chans;
    end

    % compute PC per cond 
    for condi = 1:length(conds)
        cond = conds{condi};

        % load
        disp(['Load data of condition ', upper(cond)]);
        eeg = getTrials(sub, cond, ic_struct.condParStruct, [], chans);

        % pars
        if si == 1 && condi == 1
            srate = round(1000*length(eeg.times) / (eeg.times(end) - eeg.times(1)));
            
            if stimsepOn 
                keyPnts = dsearchn(eeg.times', [ic_struct.times(1) wincut ic_struct.times(end)]');
                basePnts = keyPnts(1):keyPnts(2);
                asoPnts  = keyPnts(3):keyPnts(4);
            else 
                keyPnts = dsearchn(eeg.times', [ic_struct.times(1):winwid:ic_struct.times(end)]');
            end

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
            powerAllChans = compute_power(dataFilt, [], 0);  % baseline correction: 'z', '%', 0. 
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
        if stimsepOn
            disp('Average in baseline and after-stimulus-onset...')
            mat = zeros(size(data_ic,2), 3); % colnames: urevent, base, aso
            mat(:,1) = eeg.urevent;
            mat(:,2) = mean(data_ic(basePnts,:),1);
            mat(:,3) = mean(data_ic(asoPnts,:),1);
        else
            disp(['Average in every ', num2str(winwid), ' ms'])
            mat = zeros(size(data_ic,2), length(keyPnts)); % colnames: urevent, epoch1,..., epochn; n = length(keyPnts)-1
            mat(:,1) = eeg.urevent;
            for ei = 2:length(keyPnts)
                mat(:,ei) = mean(data_ic(keyPnts(ei-1):keyPnts(ei),:),1);
            end  % loop over epochs
        end

        % output
        varName = strcat(icaName,'_',lower(cond));
        f_path = fullfile(f_output, [icaName, ifelse(groupOn, '_overall_', '_'), num2str(icId), '_epochs', ifelse(normalizeOn, '_tNorm', '')]);
        f_temp = fullfile(f_path, [num2str(sub),'.mat']);
        disp(['Save as feature in ', f_temp])
        eval([varName,'= mat;'])
        if exist(f_path, 'dir') ~= 7
            disp('Create the path')
            mkdir(f_path)
        end           
        if exist(f_temp,'file') == 2
            save(f_temp, varName, '-append')
        else
            save(f_temp, varName)
        end
        disp('Saved')

    end  % loop over conds

end  % loop over subs

end

