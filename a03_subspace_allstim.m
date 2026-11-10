function a03_subspace_allstim(subid)

    %% set the path and create folder
    restoredefaultpath;
    addpath /seastor/a853898293/analysis/subfunction
    [~,nsensor,metapath,stimcellpath,~,behmeglog] = a01_server_setpath(subid);
     
    [~,label_scene] = a01_define_labelgroup('scene',10,behmeglog,stimcellpath,subid);
  
    % rsa folder
    c_folder = 'demean_rmtrend_subspace_allstim';
    currfile1 = dir(metapath);
    if ~ismember(c_folder,{currfile1.name})
        mkdir([metapath,'/',c_folder]);
    end
    % subject folder
    subjectpath = [metapath,'/',c_folder];
	subjectfolder = ['sub',num2str(subid)];
    currfile2 = dir(subjectpath);
    if ~ismember(subjectfolder,{currfile2.name})
        mkdir([subjectpath,'/',subjectfolder]);
    end
    
    %% load the data
	load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_goodtrialid.mat'],'cell_goodtrial');
    load([metapath,'/demean_rmtrend_dataprepare_encoding/sub',num2str(subid),'/sub_',num2str(subid),'_kpoint_31_data.mat'],'cell_raw');  % trl*sen*time
    
    timerange_whole = num2cell(1:50);
    ntime = length(timerange_whole);
    nperm = 30;
    num_pc = 2;
    
    % two pc angles
    true_cos_pcangle = nan(6,ntime,num_pc);  % g1-2 g1-3 g2-3 l1-2 l1-3 l2-3
    true_rev = nan(6,ntime);
    within_cos_pcangle = nan(6,nperm,ntime,num_pc);  % g1-1 g2-2 g3-3 l1-1 l2-2 l3-3
    within_rev = nan(6,nperm,ntime);
        
    sensor_selid = true(nsensor,1);
    for irep = 1:ntime
        tic
        timerange = timerange_whole{irep};

        pair = [1,2;1,3;2,3];
        for ipair = 1:size(pair,1)
            [true_cos_pcangle(ipair,irep,:),true_rev(ipair,irep)] = ...
                angle_stimbase_1split_true(sensor_selid,cell_raw,cell_goodtrial,'g',pair(ipair,1),pair(ipair,2),timerange,label_scene,num_pc);
            [true_cos_pcangle(3+ipair,irep,:),true_rev(3+ipair,irep)] = ...
                angle_stimbase_1split_true(sensor_selid,cell_raw,cell_goodtrial,'l',pair(ipair,1),pair(ipair,2),timerange,label_scene,num_pc);
        end

        for iperm = 1:nperm
            for igl = 1:3
                [within_cos_pcangle(igl,iperm,irep,:),within_rev(igl,iperm,irep)] = ...
                    angle_stimbase_2split_w(sensor_selid,cell_raw,cell_goodtrial,'g',igl,timerange,label_scene,num_pc);
                [within_cos_pcangle(3+igl,iperm,irep,:),within_rev(3+igl,iperm,irep)] = ...
                    angle_stimbase_2split_w(sensor_selid,cell_raw,cell_goodtrial,'l',igl,timerange,label_scene,num_pc);
            end

        end
        toc
    end
    save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_stimbase.mat'],...
        'true_cos_pcangle','true_rev','within_cos_pcangle','within_rev','-v7.3');
end