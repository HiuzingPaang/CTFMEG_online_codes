function a04_ple(subid)

     %% set the path and create folder
    restoredefaultpath;
    addpath /seastor/a853898293/analysis/subfunction
    [channel_selected,nsensor,metapath,~,subfolder,behmeglog] = a01_server_setpath(subid);  % 1 for group server
    
    % classify folder
    c_folder = 'ple_spectrum';
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
    PLE_meg_spectrum(allpre,subid,'prelearn');
    PLE_meg_spectrum(allpre,subid,'videowatching');

end