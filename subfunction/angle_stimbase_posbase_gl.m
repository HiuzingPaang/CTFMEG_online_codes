function [cos_pcangle,rev] = angle_stimbase_posbase_gl(sensor_selid,cell_raw,cell_goodtrial,label1,label2,timerange)

    num_run = size(cell_raw,1);
    colid = randperm(num_run);
    run_1 = colid(1:5); 
    run_2 = colid(6:10);

    if strcmp(label1,'g')
        wlabel_run1 = repmat([1,1,1,2,2,2,3,3,3],1,9*num_run/2)';
    elseif strcmp(label1,'l')
        wlabel_run1 = repmat([1,2,3,1,2,3,1,2,3],1,9*num_run/2)';
    end
    if strcmp(label2,'g')
        wlabel_run2 = repmat([1,1,1,2,2,2,3,3,3],1,9*num_run/2)';
    elseif strcmp(label2,'l')
        wlabel_run2 = repmat([1,2,3,1,2,3,1,2,3],1,9*num_run/2)';
    end

    run_trial = 81;
    ntime = 200;
    rawdata_run1 = nan(num_run*run_trial/2,sum(sensor_selid),ntime);
    rawdata_run2 = nan(num_run*run_trial/2,sum(sensor_selid),ntime);
    goodlabel_run1 = nan(num_run*run_trial/2,1);
    goodlabel_run2 = nan(num_run*run_trial/2,1);
    for irun = 1:num_run/2
        s = (irun-1)*run_trial + 1;
        e = irun*run_trial;
        rawdata_run1(s:e,:,:) = cell_raw{run_1(irun),1}(:,sensor_selid,:);
        rawdata_run2(s:e,:,:) = cell_raw{run_2(irun),1}(:,sensor_selid,:);
        goodlabel_run1(s:e,1) = cell_goodtrial{run_1(irun),1};
        goodlabel_run2(s:e,1) = cell_goodtrial{run_2(irun),1};
    end
    clear cell_raw cell_goodtrial

    % demean within region
    rawdata_run1 = rawdata_run1 - repmat(mean(rawdata_run1,2),1,sum(sensor_selid),1);
    rawdata_run2 = rawdata_run2 - repmat(mean(rawdata_run2,2),1,sum(sensor_selid),1);

    % data for svd
    alldata_run1 = nan(3,sum(sensor_selid),ntime);
    alldata_run2 = nan(3,sum(sensor_selid),ntime);
    for igroup = 1:3
        tmp = rawdata_run1(wlabel_run1==igroup & goodlabel_run1==1,:,:);
        alldata_run1(igroup,:,:) = mean(tmp,1);
        tmp = rawdata_run2(wlabel_run2==igroup & goodlabel_run2==1,:,:);
        alldata_run2(igroup,:,:) = mean(tmp,1);
    end
    clear tmp rawdata_run1 rawdata_run2

    num_pc = 2;
    [cos_pcangle,rev] = angle_stimbase_svd(alldata_run1,alldata_run2,timerange,num_pc);
    clear alldata_run1 alldata_run2

end