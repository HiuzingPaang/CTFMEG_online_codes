function a02_classify_enc_gl_cross_rmtest(subid,gltype,postype,num_shuffle)

% remove serial 1 5 9 trials from the test dataset
% gltype: 1=train_global_test_local; 2=train_local_test_global
% postype: 1=positions 1-3; 2=positions 2-3

%% set the path and create folder
restoredefaultpath;
addpath /seastor/a853898293/analysis/subfunction
[~,nsensor,metapath,stimcellpath,~,behmeglog] = a01_server_setpath(subid);

metapath_new = [metapath,'_revision'];
% classify folder
c_folder = 'a02_classify_enc_gl_cross_rmtest';
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

kfold = 5;
setlambda = 0.005;
rng(subid*100000+num_shuffle*10000)
each_shuffle = 1000;

%% define label_group for serial positions
numrun = size(behmeglog,2);
[~,label_seg] = a01_define_labelgroup('seg',numrun,behmeglog,stimcellpath,subid);
    
%% concatenate the data
run_trial = 81;
mat_goodtrial = nan(numrun*run_trial,1);
run_datapoint = 200;
load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_goodtrialid.mat'],'cell_goodtrial');
load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_data.mat'],'cell_raw');  % trl*sen*time
for irun = 1:numrun
    is = (irun-1)*run_trial + 1;
    ie = irun*run_trial;
    mat_goodtrial(is:ie,1) = cell_goodtrial{irun,1};
end
    
%% classification, raw signal
[cvind_group,n_comb] = a01_classify_create_groupind(mat_goodtrial,label_seg,kfold,1);
clear mat_goodtrial label_seg

whole_hit = cell(run_datapoint,1);
whole_hit_shuffle = cell(run_datapoint,1);  

for itime = 1:run_datapoint
    % extract the raw run data
    xall_raw = nan(numrun*run_trial,nsensor);
    for irun = 1:numrun
        irun_s = (irun-1)*run_trial + 1;
        irun_e = irun*run_trial;
        xall_raw(irun_s:irun_e,:) = cell_raw{irun,1}(:,:,itime);
    end
    
    % create the data of each group
    [xall_group,~] = a01_classify_create_groupdata(cvind_group,n_comb,xall_raw,nsensor);
    clear xall_raw
    
    if postype==1
        % true
        whole_hit{itime,1} = a02_classify_cvlasso_rmtest_posall(0,setlambda,n_comb,kfold,xall_group,nsensor,gltype);
        % shuffle
        shuffle_fold_hit = nan(each_shuffle,1);
        for ishuffle = 1:each_shuffle
            shuffle_fold_hit(ishuffle) = a02_classify_cvlasso_rmtest_posall(1,setlambda,n_comb,kfold,xall_group,nsensor,gltype);
        end
        whole_hit_shuffle{itime,1} = shuffle_fold_hit;
    elseif postype==2
        % true
        whole_hit{itime,1} = a02_classify_cvlasso_rmtest_pos23(0,setlambda,n_comb,kfold,xall_group,nsensor,gltype);
        % shuffle
        shuffle_fold_hit = nan(each_shuffle,1);
        for ishuffle = 1:each_shuffle
            shuffle_fold_hit(ishuffle) = a02_classify_cvlasso_rmtest_pos23(1,setlambda,n_comb,kfold,xall_group,nsensor,gltype);
        end
        whole_hit_shuffle{itime,1} = shuffle_fold_hit;
    end
end

save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),...
    '_gltype_',num2str(gltype),'_postype_',num2str(postype),'_',num2str(num_shuffle),'.mat'],...
    'whole_hit','whole_hit_shuffle','-v7.3');

end