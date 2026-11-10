function [xall_group,xall_group_trlind] = a01_classify_create_groupdata(cvind_group,n_comb,xall,nsensor)

tmp = size(cvind_group);
kfold = tmp(1);
num = tmp(2);
clear tmp

xall_group = nan(n_comb*kfold,num,nsensor);
xall_group_trlind = nan(n_comb*kfold,num);

for inum = 1:num
    for ifold = 1:kfold
        ind_trial = cvind_group{ifold,inum};
        i_s = (ifold-1)*n_comb + 1;
        i_e = ifold*n_comb;
        xall_group(i_s:i_e,inum,:) = xall(ind_trial,:);
        xall_group_trlind(i_s:i_e,inum) = ind_trial;
    end
end
clear cvind_group xall
    
end