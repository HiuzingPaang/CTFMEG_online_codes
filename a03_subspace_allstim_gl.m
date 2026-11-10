function a03_subspace_allstim_gl(subid)

    %% set the path and create folder
    restoredefaultpath;
    addpath /seastor/a853898293/analysis/subfunction
    [~,nsensor,metapath,~,~,~] = a01_server_setpath(subid);  % 1 for group server
  
    % rsa folder
    c_folder = 'demean_rmtrend_subspace_allstim_gl_3cond';
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
    
    timerange_whole = num2cell(1:200);
    ntime = length(timerange_whole);
    nperm = 30;
    
    % two pc angles
    gl_b_cos_pcangle = nan(nperm,ntime,2);  % g-l
    gl_b_rev = nan(nperm,ntime);
    gl_w_cos_pcangle = nan(2,nperm,ntime,2);  % g-g l-l
    gl_w_rev = nan(2,nperm,ntime);
        
    sensor_selid = true(nsensor,1);
    for irep = 1:ntime
        tic
        timerange = timerange_whole{irep};
        for iperm = 1:nperm
            [gl_b_cos_pcangle(iperm,irep,:),gl_b_rev(iperm,irep)] = angle_stimbase_posbase_gl(sensor_selid,cell_raw,cell_goodtrial,'g','l',timerange);
            [gl_w_cos_pcangle(1,iperm,irep,:),gl_w_rev(1,iperm,irep)] = angle_stimbase_posbase_gl(sensor_selid,cell_raw,cell_goodtrial,'g','g',timerange);
            [gl_w_cos_pcangle(2,iperm,irep,:),gl_w_rev(2,iperm,irep)] = angle_stimbase_posbase_gl(sensor_selid,cell_raw,cell_goodtrial,'l','l',timerange);
        end
        toc
    end
    save([subjectpath,'/',subjectfolder,'/sub_',num2str(subid),'_stimbase.mat'],...
        'gl_b_cos_pcangle','gl_b_rev','gl_w_cos_pcangle','gl_w_rev','-v7.3');
end