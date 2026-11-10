clear; clc; close all;

addpath D:\Matlab_Workspace\MemTim_analysis\meg_1_code_run\subfunction

[x,y] = meshgrid([1,1,1,2,2,2,3,3,3]);
xmat_g = (x==y)*1;
xmat_g(eye(9)==1) = nan;
t1 = tril(rand(9,9),-2)~=0;
t1 = t1' + t1;
xmat_g(t1==1 & (x==y)*1==1) = 0.5;

t1 = tril(rand(9,9),-3)~=0;
t1 = t1' + t1;
t2 = tril(rand(9,9),-4)~=0;
t2 = t2' + t2;
t3 = tril(rand(9,9),-6)~=0;
t3 = t3' + t3;
t4 = tril(rand(9,9),-7)~=0;
t4 = t4' + t4;
xmat_l = zeros(9,9);
xmat_l(t1==1 & t2==0) = 1;
xmat_l(t3==1 & t4==0) = 0.5;
xmat_l(eye(9)==1) = nan;

t1 = tril(rand(9,9),-1)~=0;
t1 = t1' + t1;
t2 = tril(rand(9,9),-2)~=0;
t2 = t2' + t2;
t3 = tril(rand(9,9),-3)~=0;
t3 = t3' + t3;
xmat_continuous = zeros(9,9);
xmat_continuous(t1==1 & t3==0) = 0.5;  % 0.5
xmat_continuous(t1==1 & t2==0) = 1;
xmat_continuous(eye(9)==1) = nan;

nsub = length(dir('D:\Matlab_Workspace\MemTim_behdata\expdata'))-2;
subexp = zeros(9,9,nsub);

subexp_trlwise = nan(9,180,nsub);
predict_global = nan(9,180,nsub);
predict_local = nan(9,180,nsub);
predict_serial = nan(9,180,nsub);
predict_trialtype = nan(9,180,nsub);
pos3_global = nan(26,3);
pos3_local = nan(26,3);
for isub = 1:nsub
    
    load(['D:\Matlab_Workspace\MemTim_analysis\meg_1_code_run\subfunction\WholeList_mat_v20230124\sub',...
        num2str(isub),'WholeListStimcell.mat'],'stimcell');
    ind = find(strcmp(stimcell(1,:),{'TestPos1'}));
    testpos = [cell2mat(stimcell(2:end,ind));cell2mat(stimcell(2:end,ind+1))];
    
    subdata_pre = nan(10*9,12);
    for irun = 1:10
        load(['D:\Matlab_Workspace\MemTim_behdata\expdata\sub',num2str(isub),...
            '\sub',num2str(isub),'_run',num2str(irun),'_respmat.mat'],'respmat');
        subdata_pre((irun-1)*9+1:irun*9,:) = respmat;
    end
    
    % 1=irun, 2=itrial
    % 3=rt of button1f, 4=postcur(nan,1-9,99999), 5=acc
    % 6=precur, 7=rt of moving timeline
    % 8=testpos
    subdata = [subdata_pre(:,[1,2,3,5,7,9,11]);...
        subdata_pre(:,[1,2,4,6,8,10,12])];
    % all valid trials
    crit = ~isnan(subdata(:,4)) & subdata(:,4)~=99999;
    data = [subdata(crit,:),testpos(crit,:)];
    
    mat = zeros(9,9);
    for itest = 1:9
        for ichosen = 1:9
            if itest ~= ichosen
                perc = sum(data(data(:,8)==itest,4)==ichosen)/sum(data(:,5)==0);
                mat(itest,ichosen) = perc;
            end
        end
    end
    
    subexp(:,:,isub) = mat;
    
    pos3_global_tmp = zeros(180,3);
    pos3_local_tmp = zeros(180,3);
    idx_global = [1,1,1,2,2,2,3,3,3];
    idx_local = [1,2,3,1,2,3,1,2,3];

    % for trial wise
    ic_incorr = 0;
    for iline = 1:size(data,1)
        if data(iline,4) ~= data(iline,8)  % incorrect trials
            ic_incorr = ic_incorr + 1;
            subexp_trlwise(:,ic_incorr,isub) = 0;
            subexp_trlwise(data(iline,8),ic_incorr,isub) = nan;
            subexp_trlwise(data(iline,4),ic_incorr,isub) = 1;
            predict_global(:,ic_incorr,isub) = xmat_g(data(iline,8),:);
            predict_local(:,ic_incorr,isub) = xmat_l(data(iline,8),:);
            predict_serial(:,ic_incorr,isub) = xmat_continuous(data(iline,8),:);
            predict_trialtype(:,ic_incorr,isub) = ic_incorr;  % data(iline,8)
            predict_trialtype(data(iline,8),ic_incorr,isub) = nan;

            if ismember(data(iline,4),find(xmat_g(data(iline,8),:)))  % correct global transposition error
                pos3_global_tmp(iline,idx_global(data(iline,8))) = 1;
            end
            if ismember(data(iline,4),find(xmat_l(data(iline,8),:)))  % correct global transposition error
                pos3_local_tmp(iline,idx_local(data(iline,8))) = 1;
            end
        end
    end

    pos3_global(isub,:) = sum(pos3_global_tmp,1)/sum(sum(pos3_global_tmp));
    pos3_local(isub,:) = sum(pos3_local_tmp,1)/sum(sum(pos3_local_tmp));

end

subexp_stats = nan(nsub,3);  % 1=global, 2=local, 3=continuous
subexp_stats_trlwise = nan(nsub,3);
for isub = 1:nsub
    
    y_mat = subexp(:,:,isub);
    
    % create y and predictors
    y = reshape(y_mat,[],1);
    x_global = reshape(xmat_g,[],1);
    x_local = reshape(xmat_l,[],1);
    x_seg = reshape(xmat_continuous,[],1);
    
    % remove diagonal element
    nanid = reshape(eye(9)==1,[],1);
    y = y(~nanid);
    x_global = x_global(~nanid);
    x_local = x_local(~nanid);
    x_seg = x_seg(~nanid);

    tbl = table(y,x_global,x_local,x_seg);

    m_whole = fitglm(tbl,'y ~ 1 + x_global + x_local + x_seg');
    stat.whole_r2 = m_whole.Rsquared.Ordinary;

    m_g = fitglm(tbl,'y ~ 1 + x_local + x_seg');
    stat.g_r2 = stat.whole_r2 - m_g.Rsquared.Ordinary;
    clear m_g
    
    m_l = fitglm(tbl,'y ~ 1 + x_global + x_seg');
    stat.l_r2 = stat.whole_r2 - m_l.Rsquared.Ordinary;
    clear m_l

    m_s = fitglm(tbl,'y ~ 1 + x_global + x_local');
    stat.s_r2 = stat.whole_r2 - m_s.Rsquared.Ordinary;
    clear m_s
    
    subexp_stats(isub,:) = [stat.g_r2,stat.l_r2,stat.s_r2];
    clear stat
    
    % trial wise, create y and predictors
    y = reshape(subexp_trlwise(:,:,isub),[],1);
    x_global = reshape(predict_global(:,:,isub),[],1);
    x_local = reshape(predict_local(:,:,isub),[],1);
    x_seg = reshape(predict_serial(:,:,isub),[],1);
    x_trialtype = reshape(predict_trialtype(:,:,isub),[],1);
    
    % remove nan element
    nanid = isnan(y);
    y = y(~nanid);
    x_global = x_global(~nanid);
    x_local = x_local(~nanid);
    x_seg = x_seg(~nanid);
    x_trialtype = x_trialtype(~nanid);
    
    tbl = table(y,x_global,x_local,x_seg,x_trialtype);

    m_whole = fitglme(tbl,'y ~ 1 + x_global + x_local + x_seg + (1|x_trialtype)');
    stat.whole_r2 = m_whole.Rsquared.Ordinary;

    m_g = fitglme(tbl,'y ~ 1 + x_local + x_seg + (1|x_trialtype)');
    stat.g_r2 = stat.whole_r2 - m_g.Rsquared.Ordinary;
    clear m_g
    
    m_l = fitglme(tbl,'y ~ 1 + x_global + x_seg + (1|x_trialtype)');
    stat.l_r2 = stat.whole_r2 - m_l.Rsquared.Ordinary;
    clear m_l

    m_s = fitglme(tbl,'y ~ 1 + x_global + x_local + (1|x_trialtype)');
    stat.s_r2 = stat.whole_r2 - m_s.Rsquared.Ordinary;
    clear m_s
    
    subexp_stats_trlwise(isub,:) = [stat.g_r2,stat.l_r2,stat.s_r2];
    clear stat
end

% control group trial wise

nsub = length(dir('D:\Matlab_Workspace\MemTim_behdata\condata'))-2;

subcon = nan(9,9,nsub);

subcon_trlwise = nan(9,180,nsub);
predict_global = nan(9,180,nsub);
predict_local = nan(9,180,nsub);
predict_serial = nan(9,180,nsub);
predict_topic = nan(9,180,nsub);
predict_trialtype = nan(9,180,nsub);

for isub = 1:nsub
    
    load(['D:\Matlab_Workspace\MemTim_analysis\meg_1_code_run\subfunction\WholeList_mat_v20230124_controlexp\sub',...
        num2str(isub),'WholeListStimcell_controlexp.mat'],'stimcell');
    ind = find(strcmp(stimcell(1,:),{'TestPos1'}));
    testpos = [cell2mat(stimcell(2:end,ind));cell2mat(stimcell(2:end,ind+1))];
    ind = find(strcmp(stimcell(1,:),{'LeaClip1'}));
    sceneorder = cellfun(@(x) str2double(x(34)),stimcell(2:end,ind:ind+8));
    sceneorder = repmat(sceneorder,2,1);
    
    load(['D:\Matlab_Workspace\MemTim_behdata\condata\sub',num2str(isub),...
        '\sub',num2str(isub),'respmat_controlexp.mat'],'respmat');
    subdata_pre = respmat;
    
    subdata = [subdata_pre(:,[1,2,3,5,7,9,11]);...
        subdata_pre(:,[1,2,4,6,8,10,12])];
        
    % incorrect trials
    crit = ~isnan(subdata(:,4)) & subdata(:,4)~=99999 & subdata(:,5)==0;
    data = [subdata(crit,:),testpos(crit,:)];  % 4=chosen position, 8=teseted position
    sceneorder = sceneorder(crit,:);
    
    mat = zeros(9,9);
    for itest = 1:9
        for ichosen = 1:9
            if itest ~= ichosen
                perc = sum(data(data(:,8)==itest,4)==ichosen)/sum(data(:,5)==0);
                mat(itest,ichosen) = perc;
            end
        end
    end
    subcon(:,:,isub) = mat;
    
    % for trial wise
    ic_incorr = 0;
    for iline = 1:size(data,1)
        if data(iline,4) ~= data(iline,8)  % incorrect trials
            ic_incorr = ic_incorr + 1;
            subcon_trlwise(:,ic_incorr,isub) = 0;
            subcon_trlwise(data(iline,8),ic_incorr,isub) = nan;
            subcon_trlwise(data(iline,4),ic_incorr,isub) = 1;
            predict_global(:,ic_incorr,isub) = xmat_g(data(iline,8),:);
            predict_local(:,ic_incorr,isub) = xmat_l(data(iline,8),:);
            predict_serial(:,ic_incorr,isub) = xmat_continuous(data(iline,8),:);
            
            tmp = double(sceneorder(iline,:)==sceneorder(iline,data(iline,8)));
            idx_dis = abs(find(tmp) - data(iline,8));
            [~, idx_max] = max(idx_dis);
            tmp(idx_max) = 0.5;
            predict_topic(:,ic_incorr,isub) = tmp;
            predict_topic(data(iline,8),ic_incorr,isub) = nan;
            
            predict_trialtype(:,ic_incorr,isub) = ic_incorr;  % data(iline,8)
            predict_trialtype(data(iline,8),ic_incorr,isub) = nan;
        end
    end
         
end

subcon_stats = nan(nsub,3); % 1=global, 2=local, 3=continuous
subcon_stats_trlwise = nan(nsub,4); % 4=topic
for isub = 1:nsub
    
    y_mat = subcon(:,:,isub);
    
    % create y and predictors
    y = reshape(y_mat,[],1);
    x_global = reshape(xmat_g,[],1);
    x_local = reshape(xmat_l,[],1);
    x_seg = reshape(xmat_continuous,[],1);
    
    % remove diagonal element
    nanid = reshape(eye(9)==1,[],1);
    y = y(~nanid);
    x_global = x_global(~nanid);
    x_local = x_local(~nanid);
    x_seg = x_seg(~nanid);

    tbl = table(y,x_global,x_local,x_seg);

    m_whole = fitglm(tbl,'y ~ 1 + x_global + x_local + x_seg');
    stat.whole_r2 = m_whole.Rsquared.Ordinary;

    m_g = fitglm(tbl,'y ~ 1 + x_local + x_seg');
    stat.g_r2 = stat.whole_r2 - m_g.Rsquared.Ordinary;
    clear m_g
    
    m_l = fitglm(tbl,'y ~ 1 + x_global + x_seg');
    stat.l_r2 = stat.whole_r2 - m_l.Rsquared.Ordinary;
    clear m_l

    m_s = fitglm(tbl,'y ~ 1 + x_global + x_local');
    stat.s_r2 = stat.whole_r2 - m_s.Rsquared.Ordinary;
    clear m_s
    
    subcon_stats(isub,:) = [stat.g_r2,stat.l_r2,stat.s_r2];
    clear stat
    
    % trial wise, create y and predictors
    y = reshape(subcon_trlwise(:,:,isub),[],1);
    x_global = reshape(predict_global(:,:,isub),[],1);
    x_local = reshape(predict_local(:,:,isub),[],1);
    x_seg = reshape(predict_serial(:,:,isub),[],1);
    x_topic = reshape(predict_topic(:,:,isub),[],1);
    x_trialtype = reshape(predict_trialtype(:,:,isub),[],1);
    
    % remove nan element
    nanid = isnan(y);
    y = y(~nanid);
    x_global = x_global(~nanid);
    x_local = x_local(~nanid);
    x_seg = x_seg(~nanid);
    x_topic = x_topic(~nanid);
    x_trialtype = x_trialtype(~nanid);
    
    tbl = table(y,x_global,x_local,x_seg,x_topic,x_trialtype);

    m_whole = fitglme(tbl,'y ~ 1 + x_global + x_local + x_seg + x_topic + (1|x_trialtype)');
    stat.whole_r2 = m_whole.Rsquared.Ordinary;

%     m_g = fitglm(tbl,'y ~ 1 + x_local + x_seg');
    m_g = fitglme(tbl,'y ~ 1 + x_local + x_seg + x_topic + (1|x_trialtype)');
    stat.g_r2 = stat.whole_r2 - m_g.Rsquared.Ordinary;
    clear m_g
    
    m_t = fitglme(tbl,'y ~ 1 + x_local + x_seg + x_global + (1|x_trialtype)');
    stat.t_r2 = stat.whole_r2 - m_t.Rsquared.Ordinary;
    clear m_t
    
    m_l = fitglme(tbl,'y ~ 1 + x_global + x_seg + x_topic + (1|x_trialtype)');
    stat.l_r2 = stat.whole_r2 - m_l.Rsquared.Ordinary;
    clear m_l

    m_s = fitglme(tbl,'y ~ 1 + x_global + x_local + x_topic + (1|x_trialtype)');
    stat.s_r2 = stat.whole_r2 - m_s.Rsquared.Ordinary;
    clear m_s
    
    subcon_stats_trlwise(isub,:) = [stat.g_r2,stat.l_r2,stat.s_r2,stat.t_r2];
    clear stat

end
