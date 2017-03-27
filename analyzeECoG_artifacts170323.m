load('R:\obregon\NovelvsRepeatTask\subjects\5d37ff\data\d4\5d37ff_scm\5d37ff_d4_epochs.mat')

dataR = [];
dataR.fsample       = sampfreq; % sampling rate in Hz
dataR.label         = ch_names; % each channel gets a label
dataR.time          = cell(1,length(epochs)); % each trial is represented by an element in cell arrays
dataR.trial         = cell(1,length(epochs)); % "dataR.time" & "dataR.trial" (leave empty for now)
% dataR.sampleinfo    = nan(length(epochs),2); % each row indexes the first and last sample of each trial

trllng = [];
for k=1:length(epochs)
    trllng(k) = size(epochs{k},2);
end

for k=1:length(epochs)
    
    % time values for each trial, in sec (assumes that the sampling rate is 1000 Hz)
    % each time value corresponds to a column in dataR.trial{k}
    dataR.time{k} = (0:trllng(k)-1)/sampfreq;
    
    % choose data for each trial based around event_bgn(k), which marks the trial onset
    % each column contains data for one timepoint; rows correspond to channels
    dataR.trial{k} = epochs{k}*10e6;
    
%     % write the onset and offset samples to dataR.sampleinfo
%     dataR.sampleinfo(k,:)=[event_bgn(k)-fs event_bgn(k)+(fs*10)];
    
end


c = 0;
inds = [];

for k = 1:length(epochs)
    int = c+1:c+size(epochs{k},2);
    inds(k,:) = [int(1) int(end)];
    c = c + size(epochs{k},2) + round(sampfreq);
end
   
dataR.sampleinfo = inds;

cfg = [];
cfg.demean = 'yes';
cfg.detrend = 'yes';
dataR = ft_preprocessing(cfg,dataR);

% do the spectral analysis - time-averaged
clear cfgpow
cfgpow.output      = 'pow';
cfgpow.method      = 'mtmfft';
cfgpow.pad         = 'maxperlen';
cfgpow.keeptrials  = 'no';
cfgpow.taper       = 'hanning';
cfgpow.foilim      = [1,500];
% cfgpow.foilim      = [1,100];
freqR = ft_freqanalysis(cfgpow, dataR);

cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = 60;
dataF = ft_preprocessing(cfg,dataR);

% do the spectral analysis - time-averaged
clear cfgpow
cfgpow.output      = 'pow';
cfgpow.method      = 'mtmfft';
cfgpow.pad         = 'maxperlen';
cfgpow.keeptrials  = 'no';
cfgpow.taper       = 'hanning';
cfgpow.foilim      = [1,500];
% cfgpow.foilim      = [1,100];
freqF = ft_freqanalysis(cfgpow, dataF);

figure;plot(freqR.freq,[freqR.powspctrm(1,:); freqF.powspctrm(1,:)])

cfg = [];
cfg.dftfilter = 'yes';
cfg.dftfreq = [60:60:300 59.9:60:300 60.1:60:300];
dataF2 = ft_preprocessing(cfg,dataR);

% do the spectral analysis - time-averaged
clear cfgpow
cfgpow.output      = 'pow';
cfgpow.method      = 'mtmfft';
cfgpow.pad         = 'maxperlen';
cfgpow.keeptrials  = 'no';
cfgpow.taper       = 'hanning';
cfgpow.foilim      = [1,500];
% cfgpow.foilim      = [1,100];
freqF2 = ft_freqanalysis(cfgpow, dataF2);

figure;plot(freqR.freq,[freqR.powspctrm(1,:); freqF2.powspctrm(1,:)])

% do the time-locked analysis
clear cfg
cfg.channel       = 'all';
cfg.covariance    = 'no';
cfg.keeptrials    = 'no';
cfg.removemean    = 'no';
cfg.vartrllength  =  2;
timelock=ft_timelockanalysis(cfg,dataF2);

figure;plot(timelock.time,timelock.avg)
legend(timelock.label); legend off
