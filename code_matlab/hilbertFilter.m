function dataFilt = hilbertFilter(data, srate, bandrange, checkOn)
% dataFilt = hilbertFilter(data, srate, bandrange, checkOn)
% do Hilbert-filtering 
% data: nChan x nPnt x nTrial

% set pars
nyquist = srate/2;
lower_filter_bound = min(bandrange); % Hz
upper_filter_bound = max(bandrange); % Hz
transition_width   = 0.2;  % used to distortion in time-domain
filter_order       = round(3*srate/lower_filter_bound);  % contain >= 3 cycles of the lowest freq, in unit of pnt

% build filter
ffrequencies  = [ 0 (1-transition_width)*lower_filter_bound lower_filter_bound upper_filter_bound (1+transition_width)*upper_filter_bound nyquist ]/nyquist;
idealresponse = [ 0 0 1 1 0 0 ];
filterweights = firls(filter_order, ffrequencies, idealresponse);  % create finite-impulse response filters via least squares
% the length of the returned filterweight is filter_order + 1
% fir1 is preferred in narrow-band
% filterweights = fir1(filter_order, [lower_filter_bound, upper_filter_bound]); 

% check kernel
if checkOn
    filterweightsNorm = zscore(filterweights);
    hz_filtkern = linspace(0, nyquist, floor(filter_order/2)+1);

    figure
    plot(ffrequencies*nyquist, idealresponse, 'r')
    hold on

    fft_filtkern  = abs(fft(filterweightsNorm));
    fft_filtkern  = fft_filtkern./max(fft_filtkern); % normalized to 1.0 for visual comparison ease
    plot(hz_filtkern, fft_filtkern(1:ceil(length(fft_filtkern)/2)),'b')

    set(gca,'ylim',[-.1 1.1],'xlim',[0 nyquist])
    legend({'ideal';'best fit'})

    freqsidx = dsearchn(hz_filtkern', ffrequencies'*nyquist);
    if size(fft_filtkern, 2) ~= length(fft_filtkern)
        fft_filtkern = fft_filtkern';
    end
    title(['SSE: ', num2str(sum((idealresponse - fft_filtkern(freqsidx)).^2 ))])
    % the closer SSE to zero, the better; don't use filter with SSE > 1

end

% filter data 
rawsize = size(data);

% in case data only has one trial (size(data,3) == 1)
if size(data,3) == 1 || length(rawsize) == 2
    rawsize = [rawsize 1];
end
    
if rawsize(2) <= filter_order * 3  % zero-padding if pnt number isnt enough
    data = [data zeros(rawsize(1), filter_order * 3 - rawsize(2) + 1, rawsize(3))];
end
dataFilt = zeros(size(data));
for ti = 1:size(data,3)
    for chani = 1:size(data,1)
        dataFilt(chani,:,ti) = filtfilt(filterweights, 1, double(data(chani,:,ti)));
    end 
end 
data     =     data(:,1:rawsize(2),:);
dataFilt = dataFilt(:,1:rawsize(2),:);

% transform to analytical signal
for ti = 1:size(dataFilt,3)
    dataFilt(:,:,ti) = hilbert(squeeze(dataFilt(:,:,ti))')'; % time should be in the first dimension. 
end % triali

% check transform by examing phase
if checkOn
    chan = randsample(1:rawsize(1),1);
    ni = randsample(1:rawsize(3),1);  % pick a trial randomly
    times = [0:rawsize(2)-1] * 1000/srate;
    figure
    suptitle(['Channel ', num2str(chan)])
    subplot(1,2,1)
    plot(times, data(chan,:,ni), 'b', times, real(dataFilt(chan,:,ni)), 'r--')
    xlabel('Time (ms)'), ylabel('Amp. (\muV)')
    legend({'Orignial','Filtered'})
    subplot(1,2,2)
    plot(times,angle(dataFilt(chan,:,ni)'),'b');
    xlabel('Time (ms)'), ylabel('Phase angle (rad.)')
end

end