function EEG = preprocessing(data, stage)
% preprocessing(data, stage)
% data can be multiple in some stages
% stage 1: import data, filtering, downsampling, interpolation, epoching
% stage 2: visually inspect epochs before ICA
% stage 3: remove marked epochs and run ICA
% stage 4: visually inspect components
% stage 5: remove marked comps and visually inspect epochs again
% stage 6: remove marked epochs 
% Author: Christina Jin (christina.mik109@gmail.com)

f_main = fileparts(which('mind_wandering2'));
cd(f_main)

% check input pars
if numel(stage)~=1
    error('One stage can be processed at one time!')
end

% set pars
f_chanlocs    = 'biosemi32.ced';
f_raw         = fullfile(f_main, 'raw'); % for raw data
f_pars        = 'pars_preprocessing.mat'; % for saving/loading some individual pars
srate         = 256; 
bandpass      = [0.1 42]; 
epoch         = [-1 3]; %in seconds
baseline      = [-0.2 0];
triggers      = {'11','22','33','44'}; %raw triggers; they will be re-reedited later
cmps2plot     = 32:-1:1;

% progress bar
n  = length(data);
pi = 0;
progressbar('Processing...')

for i = data
        
    if stage == 1

        % import raw
        EEG = pop_biosig(fullfile(f_raw, ['sub',num2str(i),'.bdf']), 'channels', 1:38, 'ref',[37 38],'refoptions',{'keepref' 'off'}); %import & reref; reference channels are not kept here

        % edit chan struc
        %EEG = pop_chanedit(EEG, 'lookup', f_chanlocs); %load channal info
        EEG = pop_chanedit(EEG,'load',{f_chanlocs, 'filetype' 'autodetect'});
        
        % add ref info
        EEG = pop_chanedit(EEG, 'setref', {'1:36' 'TP9 TP10'}); %set reference info: mastoids
        
        % rm extra chans
        EEG = pop_select(EEG, 'nochannel', 33:36); %remove extra channels because there is no coordinates for these channels which might affects ICA
        
        % plot to check if the location is correct
        % figure; topoplot(1:32,EEG.chanlocs,'electrodes','labelpoint','chaninfo',EEG.chaninfo);
        
        % save 
        EEG = pop_saveset(EEG,['preprocessing\\',num2str(i),'_reref.set']);

        % interpolate bad channels 
        load(f_pars, 'chancorrect')
        rowno = find([chancorrect{:,1}]==i);
        if isempty(rowno)
            bad_chan  = [];
            srnd_chan = [];
        else
            disp('Bad channels recorded')
            bad_chan  = chancorrect{rowno,2};
            srnd_chan = chancorrect{rowno,3};
        end
        for k = 1:length(bad_chan)
          bad  = bad_chan(k);
          srnd = srnd_chan(k,:);
          srnd(isnan(srnd)) = [];
          EEG  = bad_chan_correct(EEG, bad, srnd); % to correct the bad channel by interpolation
        end

        % band-pass filtering
        EEG = pop_eegfiltnew(EEG, min(bandpass), max(bandpass));

        % down-sampling
        EEG = pop_resample(EEG, srate);

        % extract epochs
        EEG = pop_epoch(EEG, triggers, epoch); 
        
        % baseline correction
        EEG = pop_rmbase(EEG, baseline*1000);
        %if EEG.trials ~= 825
        %    error('Incorrect trigger number!')
        %end
        
        % count trials
        load(f_pars, 'trialCount0_epoching')
        trialCount0_epoching(i) = EEG.trials; 
        save(f_pars, 'trialCount0_epoching','-append')

        % save 
        pop_saveset(EEG, fullfile('preprocessing',[num2str(i),'_epochs.set']));

    end

    %%
    % maybe inspect trials before ICA?
    if stage == 2

        % load data
        EEG = pop_loadset(fullfile('preprocessing',[num2str(i),'_epochs.set']));
        
        % mannual processing
        disp ('Call the following:')
        disp ('pop_eegplot(EEG);')  % visual inpection
        
        % save markers
        disp (['pop_saveset(EEG, fullfile(''preprocessing'', ''', num2str(i),'_epochs.set''))']);
        
    end

    %%

    if stage == 3

        % load 
        EEG = pop_loadset(fullfile('preprocessing', [num2str(i),'_epochs.set']));
        
        % register marked trials
        load(f_pars, 'eid2rej_beforeICA', 'uid2rej_beforeICA')
        eid2rej_beforeICA{i} = find(EEG.reject.rejmanual);
        ureventList = [EEG.event.urevent];
        uid2rej_beforeICA{i} = ureventList(eid2rej_beforeICA{i});
        save(f_pars, 'uid2rej_beforeICA', 'eid2rej_beforeICA', '-append');
        
        % reject the marked epoches
        EEG = pop_rejepoch(EEG, EEG.reject.rejmanual);
        
        % count trials
        load(f_pars, 'trialCount1_beforeICA')
        trialCount1_beforeICA(i) = EEG.trials;
        save(f_pars, 'trialCount1_beforeICA', '-append')
        
        % save before run ICA
        pop_saveset(EEG, fullfile('preprocessing', [num2str(i),'_epochs_ica.set']));

        % run ica 
        EEG = pop_runica(EEG, 'extended', 1, 'stop', 1E-7);

        % overwrite 
        pop_saveset(EEG, fullfile('preprocessing', [num2str(i),'_epochs_ica.set']));
        
    end

    %%

    if stage == 4

        % load 
        EEG = pop_loadset(fullfile('preprocessing', [num2str(i),'_epochs_ica.set']));

        % plot cmps
        disp('Call the following:')
        disp(['pop_prop(EEG, 0,[', num2str(cmps2plot),'], NaN,{''freqrange'' [0.1 42]});'])
        
        % save markers
        disp (['pop_saveset(EEG, fullfile(''preprocessing'', ''', num2str(i),'_epochs_ica.set''))']);

    end

    %%

    if stage == 5
        
        % load 
        EEG = pop_loadset(fullfile('preprocessing', [num2str(i),'_epochs_ica.set']));
     
        % register artifact comps
        load(f_pars,'rm_comp')
        rm_comp{i} = find(EEG.reject.gcompreject);
        save(f_pars,'rm_comp','-append');
        
        % remove comps
        EEG = pop_subcomp(EEG);
        
        % baseline correction
        EEG = pop_rmbase(EEG, baseline*1000);
        
        % save
        EEG = pop_saveset(EEG, fullfile('preprocessing', [num2str(i),'_epochs_ica_a.set']));
        
        % visually inspect epochs again
        disp ('Call the following:')
        disp ('pop_eegplot(EEG);')
        
        % save markers
        disp (['pop_saveset(EEG, fullfile(''preprocessing'', ''', num2str(i),'_epochs_ica_a.set''))']);

    end

    %%

    if stage == 6
        
        % load
        EEG = pop_loadset(fullfile('preprocessing', [num2str(i),'_epochs_ica_a.set']));
        
        % register marked trials
        load(f_pars, 'eid2rej_afterICA', 'uid2rej_afterICA')
        eid2rej_afterICA{i} = find(EEG.reject.rejmanual);
        ureventList = [EEG.event.urevent];
        uid2rej_afterICA{i} = ureventList(eid2rej_afterICA{i});
        save(f_pars, 'uid2rej_afterICA', 'eid2rej_afterICA', '-append');
        
        % reject the marked epoches
        EEG = pop_rejepoch(EEG, EEG.reject.rejmanual);
        
        % count trials
        load(f_pars, 'trialCount2_afterICA')
        trialCount2_afterICA(i) = EEG.trials;
        save(f_pars, 'trialCount2_afterICA', '-append')
        
        % save 
        pop_saveset(EEG, fullfile('preprocessing', [num2str(i),'_epochs_ica_a2.set']));
        
        % save urevents
        EEGraw = pop_loadset(fullfile('preprocessing',[num2str(i),'_epochs.set']));
        urevent = [EEGraw.urevent.type];
        save(fullfile('urevent_matfile', [num2str(i), '.mat']), 'urevent');

    end
    
    %%
    pi = pi + 1;
    progressbar(pi/n)
        
end % i: sub

end %func

