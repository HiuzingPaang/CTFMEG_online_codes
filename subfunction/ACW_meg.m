function ACW_meg(allpre,subid,period,window,overlap,maxlag)

subfolder = allpre.subfolder;
behmeglog = allpre.behmeglog;
nsensor = allpre.nsensor;
metapath = allpre.metapath;
channel_selected = allpre.channel_selected;
subjectpath = allpre.subjectpath;
subjectfolder = allpre.subjectfolder;

if nargin<4
    window = 20;  % second
    overlap = 0.5;  % percentage
    maxlag = 3;  % second
end

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
    sr = 1200;
    window = window*sr;
    maxlag = maxlag*sr;
    
    %% prelearn period
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
    % calculate ACW
    timescale = nan(nsensor,1);
    for isensor = 1:nsensor
        acf = [];
        iwin = 1;
        while 1
            x = mat_rawdata(:,isensor);
            t_window = [1+(iwin-1)*window*overlap,window+(iwin-1)*window*overlap];
            if t_window(2)>=length(x)
                break
            end
            % Normalizes the sequence so that the autocorrelations at zero lag equal 1
            acf(iwin,:) = xcorr(x(t_window(1):t_window(2)),maxlag,'coeff');
            iwin = iwin + 1;
        end
        acf_mean = mean(acf,1);
        [~,ind_max] = max(acf_mean);
        ind_half = acf_mean>=0.5;
        under_half = find(ind_half(ind_max:end)==0);
        first_under_half = under_half(1);
        acw = 1000*(2*(first_under_half-1)-1)/sr;  % ms
        timescale(isensor,1) = acw;
        %% save it
        save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_ACWmeg_',period,'_window_',num2str(window/sr),'.mat'],'timescale');
    end
    clear  % clear all vars

end