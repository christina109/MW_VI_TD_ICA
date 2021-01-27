function plotTrial(sub, urevent, chan, timewin)

% defaults
if nargin < 4 
    timewin = [];
end

% path
f_main = fileparts(which('mind_wandering2'));
cd(f_main)

% load
EEG = pop_loadset([f_main, 'preprocessing\\', num2str(sub), '_epochs_ica_a2.set']);

% trial id
idArr = [EEG.event.urevent] == urevent;
epochArr = [EEG.event.epoch];
tid = max(epochArr(idArr));  % in case epochs with two triggers

% chan id
idArr = strcmpi({EEG.chanlocs.labels}, chan);

% pnt id 
if isempty(timewin)
    data = EEG.data(idArr, :, tid);
else
    pidrange = dsearchn(EEG.times', timewin');
    data = EEG.data(idArr, pidrange(1):pidrange(2), tid);
end

% plot
figure
if isempty(timewin)
    plot(EEG.times, data);
else
    plot(EEG.times(pidrange(1):pidrange(2)), data);
end


end