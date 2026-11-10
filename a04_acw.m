function a04_acw(subid)

    %% set the path and create folder
    restoredefaultpath;
    addpath /seastor/a853898293/analysis/subfunction
    [channel_selected,nsensor,metapath,~,subfolder,behmeglog] = a01_server_setpath(subid);
  
    % classify folder
    c_folder = 'timescale_both20s';
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
    
    %% perpare allpre
	allpre.subfolder = subfolder;
    allpre.behmeglog = behmeglog;
    allpre.nsensor = nsensor;
    allpre.metapath = metapath;
    allpre.channel_selected = channel_selected;
    allpre.subjectpath = subjectpath;
    allpre.subjectfolder = subjectfolder;
    
    %% call the key functions
    ACW_meg(allpre,subid,'prelearn');
    ACW_meg(allpre,subid,'videowatching');
    ACW_meg(allpre,subid,'blank');
    
end