function computeBandICA(band, subs, condcell, condParStruct, idvOn, checkOn, saveFeatOn, datatype)
% computeBandICA(band, subs, condcell, condParStruct, idvOn, checkOn = 1,saveFeatOn = 0, datatype)
% Compute the specified band pass filtering and use
% power(default)/oscillation/phase data to run ICA
% individually (idvOn = 1) or overall (idvOn = 0)
% Author: Christina Jin (cyj.sci@gmail.com)

% default
if nargin < 6
    checkOn = 1;
end

if nargin < 7 || isempty(saveFeatOn) || ~ismember(saveFeatOn, [0 1])
    saveFeatOn = 0;
end

if nargin < 8 || ~ismember(datatype, {'power', 'phase', 'oscillation'})
    datatype = 'power';
end

% path
f_main = fileparts(which('mind_wandering2'));
f_output = fullfile(f_main, 'measure_matfile');
f_ica = fullfile(f_main, 'ica_matfile');
cd(f_main)

% pars
load('pars_markers.mat', 'bands')

% band definition
id = find(strcmpi({bands.name}, band));
if isempty(id)
    error('Band difinition doesn'' exist!')
end
measure = bands(id).name;
bandrange = bands(id).range;
chans = bands(id).chans;

% subset time window
baseline  = [-600 -200];
stimOn    = [0 600];

% defaults
if isempty(chans)  % default to all
    load('pars_EEG.mat', 'chanlocs')
    chans = {chanlocs(:).labels};
end

% compute PAC weigths individually
if idvOn == 1  
    
    % progress
    n = length(subs);
    i = 0;
    progressbar(['Computing ICA ', datatype, ' in ', upper(strrep(band, '_ica', '')), ' band per individual...'])
    
    % initialize 
    dt = struct();
    
    for subi = 1:length(subs)
        sub = subs(subi);

        for condi = 1:length(condcell)
            cond = condcell{condi};

            % load
            disp(['Load data of condition ', upper(cond)]);
            eeg = getTrials(sub, cond, condParStruct, [], chans);

            % combine conditions
            if condi == 1
                disp('Initilaize data matrix...')
                data = eeg.data;
            else 
                disp('Combine with previouly loaded')
                data(:,:,end+1:end+size(eeg.data,3)) = eeg.data;
            end
        end  % loop over conds

        % get pars
        if subi == 1
            times = eeg.times;
            srate = round(1000*length(eeg.times) / (eeg.times(end) - eeg.times(1)));
            baseIdx   = dsearchn(times', baseline');
            stimOnIdx = dsearchn(times', stimOn');
        end

        % filter & hilbert-transform
        disp('Hilber transform...')
        if subi == 1
            dataFilt = hilbertFilter(data, srate, bandrange, checkOn);
        else 
            dataFilt = hilbertFilter(data, srate, bandrange, 0);
        end

        % extract information
        if strcmpi(datatype, 'power')
            disp('Compute raw power...')
            powerAllChans = compute_power(dataFilt, [], 0);  % baseline correction: 'z', '%', 0. 
        elseif strcmpi(datatype, 'oscillation')
            disp('Transform to oscillation activity ...')
            powerAllChans = real(dataFilt);  % to be compatible with the previous version when there is only a 'power' option, output name keeps the same
        end

        % subset
        disp('Subset timewindow...')
        powerAllChans = powerAllChans(:,baseIdx(1):stimOnIdx(2),:);

        % ICA (runica): consult to pop_runica()
        disp('Comupte IC weights...')
        tmpdata = reshape( powerAllChans, length(chans), size(powerAllChans,2)*size(powerAllChans,3));
        tmprank = getrank(double(tmpdata(:,1:min(3000, size(tmpdata,2)))));
        tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean 
        
        if tmprank == size(tmpdata,1) 
            [icaweights, icasphere] = runica( tmpdata, 'lrate', 0.001, 'extended', 1, 'stop', 1E-7);  % using pars in preprocessing
        else 
            disp(['Data rank (' int2str(tmprank) ') is smaller than the number of channels (' int2str(size(tmpdata,1)) ').']);
            disp(['Proposed rank: ', num2str(tmprank)])
            [icaweights, icasphere] = runica( tmpdata, 'lrate', 0.001, 'pca', tmprank, 'extended', 1, 'stop', 1E-7 );
        end;
        
        icawinv = pinv(icaweights*icasphere); 
        
        % output structure
        dt.icaweights = icaweights;
        dt.icasphere  = icasphere;
        dt.icawinv    = icawinv;
        dt.data       = {};  % initialize here; later each cell stores one condition
        dt.band       = measure;
        dt.bandrange  = bandrange;
        dt.chans      = chans;
        dt.times      = times(baseIdx(1):stimOnIdx(2));
        dt.baseline   = baseline;
        dt.stimOn     = stimOn;
        dt.condcell   = condcell;
        dt.condParStruct = condParStruct;
       
        % save(1/2)
        eval([measure, '_', datatype, ' = dt;'])
        f_temp = fullfile(f_ica, [num2str(sub), '.mat']);
        disp(['Save ICA results as file ', f_temp])
        if exist(f_temp,'file') == 2
            save(f_temp, [measure, '_', datatype], '-append')
        else
            save(f_temp, [measure, '_', datatype])
        end
    
        for condi = 1:length(condcell)
            cond = condcell{condi};
            disp(['Compute ',datatype,' of condition ', upper(cond)]);

            % load
            eeg = getTrials(sub, cond, condParStruct, [], chans);

            % filter & hilbert-transform
            dataFilt = hilbertFilter(eeg.data, srate, bandrange, 0);

            % extract information
            if strcmpi(datatype, 'power')
                disp('Compute raw power...')
                powerAllChans = compute_power(dataFilt, [], 0);  % baseline correction: 'z', '%', 0. 
            elseif strcmpi(datatype, 'oscillation')
                disp('Transform to oscillation activity ...')
                powerAllChans = real(dataFilt);  % to be compatible with the previous version when there is only a 'power' option, output name keeps the same
            end

            % subset
            powerAllChans = powerAllChans(:,baseIdx(1):stimOnIdx(2),:);
            dt.data{condi} = powerAllChans; 
            
            % compute each ic act and save as features
            if saveFeatOn == 1
                for ici = 1:size(icaweights)

                    disp(['Compute COMPONENT ', num2str(ici), ' on ', upper(cond), '...'])

                    % get the loadings
                    loadings = icaweights(ici,:)*icasphere;

                    % initialize
                    mat = zeros(length(eeg.urevent), (stimOnIdx(2)-baseIdx(1)+2)); % colnames: urevent, val per pnt
                    mat(:,1) = eeg.urevent;
                    
                    % compute ica act: consult eeg_getdatact()
                    tp = loadings*powerAllChans(:,:);
                    tp = reshape(tp, size(powerAllChans,2), size(powerAllChans,3));
                    mat(:, 2:end) = tp';

                    % output
                    varName = [measure,'_',datatype, '_', lower(cond)];
                    f_feat = fullfile(f_output, [measure, '_', datatype, '_', num2str(ici)]);
                    if exist(f_feat) ~= 7
                        mkdir(f_feat);
                    end
                    f_temp = fullfile(f_feat, [num2str(sub),'.mat']);
                    disp(['Save as feature in ', f_temp])
                    eval([varName,'= mat;'])
                    if exist(f_temp,'file') == 2
                        save(f_temp, varName, '-append')
                    else
                        save(f_temp, varName)
                    end

                end  % loop over ics
            end  
            
        end  % loop over conds
        
        % save(2/2)
        eval([measure, '_', datatype, ' = dt;'])
        f_temp = fullfile(f_ica, [num2str(sub), '.mat']);
        disp(['Save ICA results as file ', f_temp])
        if exist(f_temp,'file') == 2
            save(f_temp, [measure, '_', datatype], '-append')
        else
            save(f_temp, [measure, '_', datatype])
        end

        % progress bar
        i = i+1;
        progressbar(i/n)

    end  % loop over subs  

% compute ovalall ICs
else
    disp(['Computing ICA ', datatype, ' in ', upper(strrep(band, '_ica', '')), ' band over ALL participants...'])
   
    % initialize output 
    dt = struct();  
    dt.data = {};
    
    % load and combine data
    for condi = 1:length(condcell)
        cond = condcell{condi};
        
        for subi = 1:length(subs)
            sub = subs(subi);

            % load
            disp(['Load data of condition ', upper(cond)]);
            eeg = getTrials(sub, cond, condParStruct, [], chans);
            data = eeg.data;
            
            % get pars
            if condi == 1 && subi == 1
                times = eeg.times;
                srate = round(1000*length(eeg.times) / (eeg.times(end) - eeg.times(1)));
                baseIdx   = dsearchn(times', baseline');
                stimOnIdx = dsearchn(times', stimOn');
            end

            % filter & hilbert-transform
            disp('Hilbert transform...')
            if condi == 1
                dataFilt = hilbertFilter(data, srate, bandrange, checkOn);
            else 
                dataFilt = hilbertFilter(data, srate, bandrange, 0);
            end

            % extract information
            if strcmpi(datatype, 'power')
                disp('Compute raw power...')
                powerAllChans = compute_power(dataFilt, [], 0);  % baseline correction: 'z', '%', 0. 
            elseif strcmpi(datatype, 'oscillation')
                disp('Transform to oscillation activity ...')
                powerAllChans = real(dataFilt);  % to be compatible with the previous version when there is only a 'power' option, output name keeps the same
            end

            % subset
            disp('Subset timewindow...')
            powerAllChans = powerAllChans(:,baseIdx(1):stimOnIdx(2),:);

            % save(1/2)
            disp('Save to output structure...')
            dt.data{subi, condi} = powerAllChans;

            eval([measure, '_', datatype, ' = dt;'])
            f_temp = fullfile(f_ica, ['group_', num2str(min(subs)),'_', num2str(max(subs)), '.mat']);
            disp(['Save loaded data as file ', f_temp])
            if exist(f_temp,'file') == 2
                save(f_temp, [measure, '_', datatype], '-append')
            else
                save(f_temp, [measure, '_', datatype])
            end
            
         end  % loop over subs
                    
    end  % loop over conds
    
    % combine over conditions x subs
    disp('Combine power data over conditions x subs...')
    for condi = 1:length(condcell)
        for subi = 1:length(subs)
            if condi == 1 && subi == 1
                powerAllChans = dt.data{subi, condi};
            else 
                powerAllChans(:,:, end + 1: end + size(dt.data{subi, condi},3)) = dt.data{subi, condi};
            end
            disp(['Have loaded ', num2str(size(powerAllChans,3)), ' data'])
        end
    end

    % ICA (runica): consult to pop_runica()
    disp('Comupte IC weights...')
    tmpdata = reshape( powerAllChans, length(chans), size(powerAllChans,2)*size(powerAllChans,3));
    tmprank = getrank(double(tmpdata(:,1:min(3000, size(tmpdata,2)))));
    tmpdata = tmpdata - repmat(mean(tmpdata,2), [1 size(tmpdata,2)]); % zero mean 

    if tmprank == size(tmpdata,1) 
        [icaweights, icasphere] = runica( tmpdata, 'lrate', 0.001, 'extended', 1, 'stop', 1E-7);  % using pars in preprocessing
    else 
        disp(['Data rank (' int2str(tmprank) ') is smaller than the number of channels (' int2str(size(tmpdata,1)) ').']);
        disp(['Proposed rank: ', num2str(tmprank)])
        [icaweights, icasphere] = runica( tmpdata, 'lrate', 0.001, 'pca', tmprank, 'extended', 1, 'stop', 1E-7 );
    end;

    icawinv = pinv(icaweights*icasphere); 

    % output structure
    dt.icaweights = icaweights;
    dt.icasphere  = icasphere;
    dt.icawinv    = icawinv;
    dt.band       = measure;
    dt.bandrange  = bandrange;
    dt.chans      = chans;
    dt.times      = times(baseIdx(1):stimOnIdx(2));
    dt.baseline   = baseline;
    dt.stimOn     = stimOn;
    dt.condcell   = condcell;
    dt.condParStruct = condParStruct;
    dt.subs       = subs;

    % save(2/2)
    eval([measure, '_', datatype, ' = dt;'])
    f_temp = fullfile(f_ica, ['group_', num2str(min(subs)),'_', num2str(max(subs)), '.mat']);
    disp(['Save ICA results as file ', f_temp])
    if exist(f_temp,'file') == 2
        save(f_temp, [measure, '_', datatype], '-append')
    else
        save(f_temp, [measure, '_', datatype])
    end
    
    % compute each IC and save as features
    if saveFeatOn 
        for subi = 1:length(subs)
            sub = subs(subi);

            for condi = 1:length(condcell)
                cond = condcell{condi};
                disp(['Compute ', datatype, ' of condition ', upper(cond)]);

                % load
                eeg = getTrials(sub, cond, condParStruct, [], chans);

                % filter & hilbert-transform
                dataFilt = hilbertFilter(eeg.data, srate, bandrange, 0);

                % extract information
                if strcmpi(datatype, 'power')
                    disp('Compute raw power...')
                    powerAllChans = compute_power(dataFilt, [], 0);  % baseline correction: 'z', '%', 0. 
                elseif strcmpi(datatype, 'oscillation')
                    disp('Transform to oscillation activity ...')
                    powerAllChans = real(dataFilt);  % to be compatible with the previous version when there is only a 'power' option, output name keeps the same
                end

                % subset
                powerAllChans = powerAllChans(:,baseIdx(1):stimOnIdx(2),:);

                for ici = 1:size(icaweights,1)

                    disp(['Compute COMPONENT ', num2str(ici), ' on ', upper(cond), '...'])
                    
                    % get the loadings
                    loadings = icaweights(ici,:)*icasphere;

                    % initialize
                    mat = zeros(length(eeg.urevent), (stimOnIdx(2)-baseIdx(1)+2)); % colnames: urevent, val per pnt
                    mat(:,1) = eeg.urevent;
                    
                    % compute ica act: consult eeg_getdatact()
                    tp = loadings*powerAllChans(:,:);
                    tp = reshape(tp, size(powerAllChans,2), size(powerAllChans,3));
                    mat(:, 2:end) = tp';

                    % output
                    varName = [measure, '_', datatype];
                    f_feat = fullfile(f_output, [measure, '_', datatype, '_overall_', num2str(ici)]);
                    if exist(f_feat, 'dir') ~= 7
                        mkdir(f_feat);
                    end
                    f_temp = fullfile(f_feat, [num2str(sub),'.mat']);
                    disp(['Save as feature in ', f_temp])
                    eval([varName,'= mat;'])
                    if exist(f_temp,'file') == 2
                        save(f_temp, varName, '-append')
                    else
                        save(f_temp, varName)
                    end

                end  % loop over pcs

            end  % loop over conds
        end  % loop over subs
    end

end 

end  % main func



function tmprank2 = getrank(tmpdata)
    % copy from pop_runica()
    
    tmprank = rank(tmpdata);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Here: alternate computation of the rank by Sven Hoffman
    %tmprank = rank(tmpdata(:,1:min(3000, size(tmpdata,2)))); old code
    covarianceMatrix = cov(tmpdata', 1);
    [~, D] = eig (covarianceMatrix);
    rankTolerance = 1e-7;
    tmprank2=sum (diag (D) > rankTolerance);
    if tmprank ~= tmprank2
        fprintf('Warning: fixing rank computation inconsistency (%d vs %d) most likely because running under Linux 64-bit Matlab\n', tmprank, tmprank2);
        tmprank2 = max(tmprank, tmprank2);
    end;
end   % getrank()