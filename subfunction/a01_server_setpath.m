function [channel_selected,nsensor,metapath,stimcellpath,subfolder,behmeglog] = a01_server_setpath(subid)

nsensor = 272;

addpath(genpath('/seastor/a853898293/addtoolboxes/eeglab2021.1/functions'));
addpath /seastor/a853898293/addtoolboxes/fieldtrip-20220113
ft_defaults;  % it takes several minutes

load('/seastor/a853898293/analysis/subfunction/cas_ctf_meglabel.mat','cas_ctf_meglabel');
channel_selected = cas_ctf_meglabel;  % exclude three bad sensors

sublist = '/seastor/a853898293/meg_cas';
metapath = '/seastor/a853898293/meg_cas_results';
stimcellpath = '/seastor/a853898293/analysis/subfunction/WholeList_mat_v20230124';

subfolder = dir([sublist,'/sub',num2str(subid),'_*']);
load([subfolder.folder,filesep,subfolder.name,filesep,'behmeglog.mat'],'behmeglog');
    
end