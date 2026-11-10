function [cos_pcangle,rev] = angle_stimbase_2split_w(sensor_selid,cell_raw,cell_goodtrial,label,labelnum,timerange,label_scene,num_pc)
    
    num_run = 10;
    if strcmp(label,'g')
        group_label = repmat([1,1,1,2,2,2,3,3,3],1,9*num_run)';
        group_o_label = repmat([1,2,3,1,2,3,1,2,3],1,9*num_run)';
    elseif strcmp(label,'l')
        group_label = repmat([1,2,3,1,2,3,1,2,3],1,9*num_run)';
        group_o_label = repmat([1,1,1,2,2,2,3,3,3],1,9*num_run)';
    end
    
    run_trial = 81;
    ntime = 200;
    rawdata_run = nan(num_run*run_trial,sum(sensor_selid),ntime);
    goodlabel_run = nan(num_run*run_trial,1);
    for irun = 1:num_run
        s = (irun-1)*run_trial + 1;
        e = irun*run_trial;
        rawdata_run(s:e,:,:) = cell_raw{irun,1}(:,sensor_selid,:);
        goodlabel_run(s:e,1) = cell_goodtrial{irun,1};
    end
    clear cell_raw cell_goodtrial

    % demean within region
    rawdata_run = rawdata_run - repmat(mean(rawdata_run,2),1,sum(sensor_selid),1);

    % data for svd
    alldata_run1 = nan(9,sum(sensor_selid),ntime);
    alldata_run2 = nan(9,sum(sensor_selid),ntime);
    for iscene = 1:3
        for ipos = 1:3
            tmp = rawdata_run(group_label==labelnum & group_o_label==ipos & label_scene==iscene & goodlabel_run==1,:,:);
            num_seltrl = size(tmp,1);
            alldata_run1((iscene-1)*3+ipos,:,:) = mean(tmp(1:floor(num_seltrl/2),:,:),1);
            alldata_run2((iscene-1)*3+ipos,:,:) = mean(tmp(floor(num_seltrl/2)+1:end,:,:),1);
        end
    end
    clear tmp rawdata_run

    [cos_pcangle,rev] = angle_stimbase_svd(alldata_run1,alldata_run2,timerange,num_pc);
    clear alldata_run1 alldata_run2
end