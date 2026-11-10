function a04_rsa_enc_4cond_gl(subid,igl,imulti)

% igl: 1=global, 2=local
% imulti: 1-7 = [10,15,20,25,30,35,40]

    %% set the path and create folder
    restoredefaultpath;
    addpath /seastor/a853898293/analysis/subfunction
    [~,~,metapath,stimcellpath,~,behmeglog] = a01_server_setpath(subid);
    
    load('/seastor/a853898293/analysis/gl_sigsensor_eachsub_oldclassify_loo.mat','subtopo_multinum');

    metapath_new = [metapath,'_revision'];
    % classify folder
    c_folder = 'a04_rsa_enc_4cond_gl';
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
    [num_seg,label_seg] = a01_define_labelgroup('seg',numrun,behmeglog,stimcellpath,subid);
    [num_g,label_g] = a01_define_labelgroup('global',numrun,behmeglog,stimcellpath,subid);
    [num_l,label_l] = a01_define_labelgroup('local',numrun,behmeglog,stimcellpath,subid);
    
    %% concatenate the data
    run_trial = 81;
    runpair = nchoosek(1:numrun,2);
    run_datapoint = 200;
    
	%% load data
    load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_goodtrialid.mat'],'cell_goodtrial');
    load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_data.mat'],'cell_raw');  % trl*sen*time

    %% rsa, sth specific information, raw power
    whole_rho_w_b = nan(run_datapoint,4);  % w_seg w_global w_local b
    sensor_sel = logical(subtopo_multinum(:,subid,igl,imulti));
    for itime = 1:run_datapoint
        A = cell(size(runpair,1),4);  % 1=within seg, 2=within global, 3=within local, 4=between
        for ipair = 1:size(runpair,1)
            run1 = cell_raw{runpair(ipair,1),1}(:,sensor_sel,itime);  % trl*sensor
            run2 = cell_raw{runpair(ipair,2),1}(:,sensor_sel,itime);
            sim_run12 = 1-pdist2(run2,run1,'spearman');
            
            run1_label_good = cell_goodtrial{runpair(ipair,1),1};
            run2_label_good = cell_goodtrial{runpair(ipair,2),1};
            sim_run12(:,run1_label_good~=1) = nan;
            sim_run12(run2_label_good~=1,:) = nan;
                
            is1 = (runpair(ipair,1)-1)*run_trial + 1;
            ie1 = runpair(ipair,1)*run_trial;
            is2 = (runpair(ipair,2)-1)*run_trial + 1;
            ie2 = runpair(ipair,2)*run_trial;
            run1_label_group_seg = label_seg(is1:ie1,1);
            run2_label_group_seg = label_seg(is2:ie2,1);
                
            run1_label_group_g = label_g(is1:ie1,1);
            run2_label_group_g = label_g(is2:ie2,1);
                
            run1_label_group_l = label_l(is1:ie1,1);
            run2_label_group_l = label_l(is2:ie2,1);
                
            W = [];
            for igroup = 1:num_seg
                w_ele = sim_run12(run2_label_group_seg==igroup,run1_label_group_seg==igroup);
                if numel(w_ele)>0
                    W = [W;reshape(w_ele,[],1)];
                    sim_run12(run2_label_group_seg==igroup,run1_label_group_seg==igroup) = nan;
                end
            end
            A{ipair,1} = W;
                
            W = [];
            for igroup = 1:num_g
                w_ele = sim_run12(run2_label_group_g==igroup,run1_label_group_g==igroup);
                if numel(w_ele)>0
                    W = [W;reshape(w_ele,[],1)];
                    sim_run12(run2_label_group_g==igroup,run1_label_group_g==igroup) = nan;
                end
            end
            A{ipair,2} = W;
                
            W = [];
            for igroup = 1:num_l
                w_ele = sim_run12(run2_label_group_l==igroup,run1_label_group_l==igroup);
                if numel(w_ele)>0
                    W = [W;reshape(w_ele,[],1)];
                    sim_run12(run2_label_group_l==igroup,run1_label_group_l==igroup) = nan;
                end
            end
            A{ipair,3} = W;
            
            B = reshape(sim_run12,[],1);
            
            A{ipair,4} = B;
        end
            
        % average them
        w_con_seg = [];
        w_con_g = [];
        w_con_l = [];
        b_con = [];
        for i = 1:size(A,1)
            w_con_seg = [w_con_seg;A{i,1}];
            w_con_g = [w_con_g;A{i,2}];
            w_con_l = [w_con_l;A{i,3}];
            b_con = [b_con;A{i,4}];
        end
        whole_rho_w_b(itime,1) = nanmean(w_con_seg);
        whole_rho_w_b(itime,2) = nanmean(w_con_g);
        whole_rho_w_b(itime,3) = nanmean(w_con_l);
        whole_rho_w_b(itime,4) = nanmean(b_con);
            
    end
    %% save it
    save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_gl_',num2str(igl),'_multinum_',num2str(imulti),'.mat'],'whole_rho_w_b');
    
end