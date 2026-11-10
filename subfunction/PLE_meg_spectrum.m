function PLE_meg_spectrum(allpre,subid,period,window)

subfolder = allpre.subfolder;
behmeglog = allpre.behmeglog;
nsensor = allpre.nsensor;
metapath = allpre.metapath;
channel_selected = allpre.channel_selected;
subjectpath = allpre.subjectpath;
subjectfolder = allpre.subjectfolder;

if nargin<4
    window = 20;  % second
end
overlap = [];  % [] is equal of 50% by default
nfft = 120000;  % points of PSD (the higher the nttf, the higher the frequency resolution, but higher noise/variability)
sr = 1200;
cutfreq = 100;

if strcmp(period,'prelearn')
    EventValue = 154;
    PreStim = 0;
    PostStim = 0.4;
    run_trial = 9;
elseif strcmp(period,'videowatching')
    EventValue = [2,4];
    PreStim = 0;
    PostStim = 2;
    run_trial = 81;
elseif strcmp(period,'blank')
    EventValue = 8;
    PreStim = 0;
    PostStim = 6;
    run_trial = 9;
end
time_pertrial = max(PreStim,PostStim);

    %% all periods
    numrun = size(behmeglog,2);
    
    mat_goodtrial = nan(numrun*run_trial,1);
    mat_rawdata = nan(time_pertrial*sr,numrun*run_trial,nsensor);
	for irun = 1:numrun
        %% obtain trl matrix
        % text_learn = 154
        % clip onset = 2
        % clip seg = 4
        % blank onset = 8
        if subid == 21
            submegfolder = dir([subfolder.folder,filesep,subfolder.name,'/meg/*',num2str(behmeglog(2,irun)),'.ds']);
        else
            submegfolder = dir([subfolder.folder,filesep,subfolder.name,'/meg/*run',num2str(behmeglog(2,irun)),'.ds']);
        end
        cfg = [];
        cfg.dataset = [submegfolder.folder,filesep,submegfolder.name];
        pre_dataname = cfg.dataset;
        cfg.trialdef.eventtype = 'UPPT001';
        cfg.trialdef.eventvalue = EventValue;
        cfg.trialdef.prestim = PreStim;  % latency in seconds
        cfg.trialdef.poststim = PostStim;  % latency in seconds
        cfg = ft_definetrial(cfg);
        trl_all = cfg.trl(:,1:3);
        clear cfg;
        
        %% define trial and cut
        data_semiclean = meg_pre_semiclean(subid,irun,pre_dataname,metapath,channel_selected);  % it only has 272 MEG channels
        cfg = [];
        cfg.trl = trl_all;
        data_semiclean_trl = ft_redefinetrial(cfg,data_semiclean);
        clear data_semiclean trl_all;
        
        %% concatenate the good/bad trial index
        artifolder = 'meg_artimark';
        load([metapath,'/',artifolder,'/sub',num2str(subid),'/jump_sub_',num2str(subid),'_run_',num2str(irun),'.mat'],'artifact_jump');
        load([metapath,'/',artifolder,'/sub',num2str(subid),'/muscle_sub_',num2str(subid),'_run_',num2str(irun),'.mat'],'artifact_muscle');
        artifact_whole = [artifact_jump;artifact_muscle];
        goodtrialid = ones(run_trial,1);
        for itrial = 1:run_trial
            for iarti = 1:size(artifact_whole,1)
                if max(data_semiclean_trl.sampleinfo(itrial,1),artifact_whole(iarti,1)) ...
                        < min(data_semiclean_trl.sampleinfo(itrial,2),artifact_whole(iarti,2))
                    goodtrialid(itrial,1) = 0;
                    break
                end
            end
        end
        irun_s = (irun-1)*run_trial + 1;
        irun_e = irun*run_trial;
        mat_goodtrial(irun_s:irun_e,1) = goodtrialid;
        
        %% concatenate the trial data
        raw_trl = nan(time_pertrial*sr,run_trial,nsensor);
        for itrial = 1:run_trial
            raw_trl(:,itrial,:) = data_semiclean_trl.trial{1,itrial}';  % chan*datapoint --> datapoint*chan
        end
        i_s = (irun-1)*run_trial + 1;
        i_e = irun*run_trial;
        mat_rawdata(:,i_s:i_e,:) = raw_trl;
        clear data_semiclean_trl raw_trl;
	end
    % good trials
    mat_rawdata = mat_rawdata(:,mat_goodtrial==1,:);
    % concatenate all good trials, AllDatapoint*chan
    mat_rawdata = squeeze(reshape(mat_rawdata,(time_pertrial*sr)*size(mat_rawdata,2),1,nsensor));
    
    % calculate power law exponent
    ple_norm = nan(nsensor,1);
    max_freq = sr/2;  % e.g., 600Hz
    step_freq = sr/nfft;  % e.g., 1200/120000 = 0.01Hz
    range_freq = 0:step_freq:max_freq;
    num_resize = 5;
    ple_norm_spectrum = nan(nsensor,sum(range_freq>0 & range_freq<=cutfreq)/num_resize);
    ple_norm_resize = nan(nsensor,1);
    
    % psd using pwelch
    [psd_series,f] = pwelch(mat_rawdata,hamming(window*sr),overlap,nfft,sr);
    psd_series_norm = psd_series./sum(psd_series,1);
    clear psd_series
    
    psd_series_norm = psd_series_norm(f<=cutfreq & f>0,:);
    f = f(f<=cutfreq & f>0);
    
    for isensor = 1:nsensor
        func = polyfit(log10(f),log10(psd_series_norm(:,isensor)),1);
        ple_norm(isensor,1) = -1*func(1);
        
        ple_norm_spectrum(isensor,:) = psd_series_norm(1:num_resize:end,isensor);
        func = polyfit(log10(f(1:num_resize:end)),log10(psd_series_norm(1:num_resize:end,isensor)),1);
        ple_norm_resize(isensor,1) = -1*func(1);
    end
    
    save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_period_',period,'_window_',num2str(window),'_ple.mat'],...
        'ple_norm','ple_norm_resize','ple_norm_spectrum','-v7.3');
    
end