function a05_boundary_subspace(subid,imulti)

% imulti: 1-7 = [10,15,20,25,30,35,40]

%% set the path and create folder
restoredefaultpath;
addpath /seastor/a853898293/analysis/subfunction
[~,nsensor,metapath,~,~,behmeglog] = a01_server_setpath(subid);

load('/seastor/a853898293/analysis/gl_sigsensor_eachsub_oldclassify_loo.mat','subtopo_multinum');
sensor_sel = logical(squeeze(subtopo_multinum(:,subid,:,imulti)));

metapath_new = [metapath,'_revision'];
% classify folder
c_folder = 'a05_boundary_subspace';
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
    
%% tidy data
numrun = size(behmeglog,2);
label_seg = repmat(1:9,1,9)';
label_no = (1:9)';
label_soft = [1,2;2,3;4,5;5,6;7,8;8,9];  % 1=pre boundary,2=post boundary
label_hard = [3,4;6,7];
    
load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_goodtrialid.mat'],'cell_goodtrial');
load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_data.mat'],'cell_raw');  % trl*sen*time
    
pre_post_time_no = 61:180;
pre_time_soft_hard = 161:200;
post_time_soft_hard = 1:80;
    
data_alltype = nan(3,nsensor,length(pre_post_time_no));  % 1=no,2=soft,3=hard
for itype = 1:3  % 1=no,2=soft,3=hard
    data = cell(numrun,1);
    if itype==1
        for irun = 1:numrun
            rawdata = cell_raw{irun,1};
            goodid = cell_goodtrial{irun,1};

            tmp = cell(size(label_no,1),1);
            for ilabel = 1:size(label_no,1)
                tmp{ilabel,1} = rawdata(label_seg==label_no(ilabel) & goodid==1,:,pre_post_time_no);
            end
            cat_str = concatenate_string(size(label_no,1),'tmp');
            tmp = eval(cat_str);  % trl*sensor*time
            data{irun,1} = tmp;
            clear tmp
        end
    elseif itype==2
        for irun = 1:numrun
            rawdata = cell_raw{irun,1};
            goodid = cell_goodtrial{irun,1};
            
            tmp = cell(size(label_soft,1),1);
            for ilabel = 1:size(label_soft,1)
                tmp_1_data = rawdata(label_seg==label_soft(ilabel,1),:,pre_time_soft_hard);
                tmp_2_data = rawdata(label_seg==label_soft(ilabel,2),:,post_time_soft_hard);
                tmp_12_data = cat(3,tmp_1_data,tmp_2_data);  % trl*sensor*time
                tmp_12_good = sum([goodid(label_seg==label_soft(ilabel,1)),goodid(label_seg==label_soft(ilabel,2))],2)==2;
                if sum(tmp_12_good)>0
                    tmp{ilabel,1} = tmp_12_data(tmp_12_good,:,:);
                end
                clear tmp_1_data tmp_2_data tmp_12_data tmp_12_good
            end
            cat_str = concatenate_string(size(label_soft,1),'tmp');
            tmp = eval(cat_str);  % trl*sensor*time
            data{irun,1} = tmp;
            clear tmp
        end
    elseif itype==3
        for irun = 1:numrun
            rawdata = cell_raw{irun,1};
            goodid = cell_goodtrial{irun,1};
            
            tmp = cell(size(label_hard,1),1);
            for ilabel = 1:size(label_hard,1)
                tmp_1_data = rawdata(label_seg==label_hard(ilabel,1),:,pre_time_soft_hard);
                tmp_2_data = rawdata(label_seg==label_hard(ilabel,2),:,post_time_soft_hard);
                tmp_12_data = cat(3,tmp_1_data,tmp_2_data);  % trl*sensor*time
                tmp_12_good = sum([goodid(label_seg==label_hard(ilabel,1)),goodid(label_seg==label_hard(ilabel,2))],2)==2;
                if sum(tmp_12_good)>0
                    tmp{ilabel,1} = tmp_12_data(tmp_12_good,:,:);
                end
                clear tmp_1_data tmp_2_data tmp_12_data tmp_12_good
            end
            cat_str = concatenate_string(size(label_hard,1),'tmp');
            tmp = eval(cat_str);  % trl*sensor*time
            data{irun,1} = tmp;
            clear tmp
        end
    end
        
    cat_str = concatenate_string(numrun,'data');
    data = eval(cat_str);  % trl*sensor*time
    data_alltype(itype,:,:) = mean(data,1);
    clear data
end
clear cell_raw cell_goodtrial
    
%% svd
num_pc = 2;
glname = {'g','l'};
for igl = 1:length(glname)
    data_sel = data_alltype(:,sensor_sel(:,igl),:);  % condition*sensor_sel*time
        
    % averaged across the pre-boundary time, for subspace estimation
    pre_mean = mean(data_sel(:,:,1:length(pre_time_soft_hard)),3);
    pre_mean = pre_mean - mean(pre_mean,1);  % condition*sensor_sel
    [U,~,~] = svd((pre_mean'*pre_mean)/(size(pre_mean,1)-1));
        
    % project
    data_pc_time = nan(size(data_sel,1),num_pc,size(data_sel,3));  % condition*pc*time
    for itime = 1:size(data_sel,3)
        data_sel_demean = data_sel(:,:,itime) - mean(data_sel(:,:,itime),1);
        data_pc_time(:,:,itime) = data_sel_demean*U(:,1:num_pc);
    end
        
    % use the boundary time point as ref
    bench = data_pc_time(:,:,length(pre_time_soft_hard)+post_time_soft_hard(1));
        
    mdd_pc = nan(size(data_pc_time,1),size(data_pc_time,3));
    for icond = 1:size(data_pc_time,1)
        tmp = squeeze(data_pc_time(icond,:,:));  % sensor_sel*time
        mdd_pc(icond,:) = sqrt(sum((tmp - bench(icond,:)').^2,1));  % euclidean distance 
    end
        
    % save
    save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_',glname{igl},'_multi_',num2str(imulti),'.mat'],'data_pc_time','mdd_pc');
end

end

function cat_str = concatenate_string(str_num,data_name)
% example: data_name -> 'w_serial_cell'
    for iserial = 1:str_num
        if iserial == 1
            cat_str = ['cat(1,',data_name,'{1}'];
        else
            cat_str = [cat_str,',',data_name,'{',num2str(iserial),'}'];
        end    
    end
    cat_str = [cat_str,')'];
end