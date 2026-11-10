function a06_content(subid,isensor)

% at each time point
% reorder 810 trials
% one subject ~ the other subjects (averaged)

    %% set the path and create folder
    restoredefaultpath;
    addpath /seastor/a853898293/analysis/subfunction
    [~,~,metapath,stimcellpath,~,behmeglog] = a01_server_setpath(subid);
    load('/seastor/a853898293/analysis/glsensor_top25.mat','glsensor');
    load('/seastor/a853898293/analysis/sub_reorder_trl.mat','sub_reorder_trl');
    
    metapath_new = [metapath,'_revision'];
    % classify folder
    c_folder = 'a06_content';
    currfile1 = dir(metapath_new);
    if ~ismember(c_folder,{currfile1.name})
        mkdir([metapath_new,'/',c_folder]);
    end
    % subject folder
    subjectpath = [metapath_new,'/',c_folder];
    subjectfolder = ['sub',num2str(subid)];
    currfile2 = dir(subjectpath);
    if ~ismember(subjectfolder,{currfile2.name})
        mkdir([subjectpath,'/',subjectfolder]);
    end

    %% define label_group
    numrun = size(behmeglog,2);
    [~,label_scene] = a01_define_labelgroup('scene',numrun,behmeglog,stimcellpath,subid);
    [~,label_animal] = a01_define_labelgroup('animal',numrun,behmeglog,stimcellpath,subid);
    [~,reid_self] = sort(sub_reorder_trl(:,subid));
    label_scene = label_scene(reid_self,1);
    label_animal = label_animal(reid_self,1);
    
    %% concatenate the data
    run_trial = 81;
    run_datapoint = 200;
    
	%% load data and compute
    data_self = nan(run_trial*numrun,sum(glsensor(:,isensor)),run_datapoint);
    load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_goodtrialid.mat'],'cell_goodtrial');
    load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_data.mat'],'cell_raw');  % trl*sen*time
    mat_goodtrl = nan(run_trial*numrun,1);
    for irun = 1:numrun
        s = (irun-1)*run_trial+1;
        e = irun*run_trial;
        mat_goodtrl(s:e,1) = cell_goodtrial{irun,1};
        data_self(s:e,:,:) = cell_raw{irun,1}(:,logical(glsensor(:,isensor)),:);
    end
    clear cell_goodtrial cell_raw
    mat_goodtrl = mat_goodtrl(reid_self,1);
    data_self = data_self(reid_self,:,:);
    
    data_other = nan(run_trial*numrun,sum(glsensor(:,isensor)),run_datapoint,25);
    sub_other = 1:26;
    sub_other(subid) = [];
    for isub_other = 1:length(sub_other)
        load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(sub_other(isub_other)),'/sub_',num2str(sub_other(isub_other)),...
            '_kpoint_31_data.mat'],'cell_raw');  % trl*sen*time
        load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(sub_other(isub_other)),'/sub_',num2str(sub_other(isub_other)),...
            '_kpoint_31_goodtrialid.mat'],'cell_goodtrial');
        tmp = nan(run_trial*numrun,sum(glsensor(:,isensor)),run_datapoint);
        tmp_goodtrl = nan(run_trial*numrun,1);
        for irun = 1:numrun
            s = (irun-1)*run_trial+1;
            e = irun*run_trial;
            tmp(s:e,:,:) = cell_raw{irun,1}(:,logical(glsensor(:,isensor)),:);
            tmp_goodtrl(s:e,1) = cell_goodtrial{irun,1};
        end
        tmp(tmp_goodtrl~=1,:,:) = nan;
        [~,reid_other] = sort(sub_reorder_trl(:,sub_other(isub_other)));
        data_other(:,:,:,isub_other) = tmp(reid_other,:,:);
        
        clear tmp cell_raw cell_goodtrial
    end
    
    data_other = nanmean(data_other,4);
    
    sim_4cond_all = nan(run_datapoint,run_datapoint,4);  % 1=video 2=animal 3=scene 4=other
    for itime1 = 1:run_datapoint
        tic
        for itime2 = 1:run_datapoint
            r_whole = 1 - pdist2(data_self(:,:,itime1),data_other(:,:,itime2),'spearman');
            r_whole(mat_goodtrl~=1,:) = nan;  % remove bad trials
        
            % video level
            sim_4cond_all(itime1,itime2,1) = nanmean(r_whole(eye(run_trial*numrun)==1));
            r_whole(eye(run_trial*numrun)==1) = nan;
            
            % animal level
            [x,y] = meshgrid(label_animal,label_animal);
            sim_4cond_all(itime1,itime2,2) = nanmean(r_whole(x==y));
            r_whole(x==y) = nan;
            
            % scene level
            [x,y] = meshgrid(label_scene,label_scene);
            sim_4cond_all(itime1,itime2,3) = nanmean(r_whole(x==y));
            r_whole(x==y) = nan;
            
            % other
            sim_4cond_all(itime1,itime2,4) = nanmean(r_whole(:));
        end
        toc
    end

    save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_event_4cond_gl_',num2str(isensor),'.mat'],'sim_4cond_all','-v7.3');
    
end