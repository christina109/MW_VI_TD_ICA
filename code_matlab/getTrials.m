function eeg = getTrials(sub, cond, condParsStruc, timewin, chans)
% eeg = getTrials(sub, cond, pars_cond, timewin = [all], chans = [all])
% timewin: in ms
% eeg.data: nChan x nPnt x nTrial

f_main = fileparts(which('mind_wandering2'));
cd(f_main)

% load pars
load('pars_EEG', 'times', 'chanlocs')

% chans id
if nargin < 5 || isempty(chans)
    chans = {chanlocs.labels};
    chansid = 1:length(chanlocs);
else
    [~,chansid] = ismember(lower(chans), lower({chanlocs.labels}));
end

% subset times
if nargin < 4 || isempty(timewin)
    subtimes = times;
    timewinid = 1:length(times);
else 
    timewinid = dsearchn(times', timewin');
    subtimes = times(timewinid(1):timewinid(end));
end

% cond id
condId = strcmpi(cond,{condParsStruc.cond});
behLabOn = condParsStruc(condId).behLabOn;
if behLabOn
    class = condParsStruc(condId).class;
    probeOn = condParsStruc(condId).probeOn;
    if probeOn 
        nback = 3;
        disp(['Load ', num2str(nback), ' trials before probes'])
    end
end

% markers
if ~behLabOn
    markers = condParsStruc(condId).marker;
else
    load(fullfile('beh_matfile', [num2str(sub), '.mat']), 'urevent', class, 'dis2pr');
    eval(['class = ', class, ';'])
end

% load
EEG = pop_loadset(['preprocessing\\', num2str(sub), '_epochs_ica_a2.set']);
    
% event id in EEG.event
[~, eventsallid] = unique([EEG.event.epoch]); % suppose first event in each epoch is time = 0

% check true event count
if length(eventsallid) ~= EEG.trials
    error('Incorrect event structure!')
end

% filter for true trigger seq
eventsall = [EEG.event.type];
eventsall = eventsall(eventsallid);

% get trial id
if ~behLabOn 
    trialid = ismember(eventsall, markers);
else
    
    if probeOn
        ur2include = urevent(ismember(class, cond) & dis2pr < nback);
    else
        ur2include = urevent(ismember(class, cond));
    end
        
    trialid = ismember([EEG.event(eventsallid).urevent], ur2include);
end

% subset
eeg.data = EEG.data(chansid, min(timewinid):max(timewinid), trialid);
eeg.urevent = [EEG.event(eventsallid(trialid)).urevent];

% add info
eeg.sub   = sub;
eeg.times = subtimes;   
eeg.chans = chans;
eeg.cond  = cond;
if behLabOn && probeOn 
    eeg.nback = nback;
end


end