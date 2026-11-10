function [cvind_group,n_comb] = a01_classify_create_groupind(label_good,label_seg,kfold,trltype)
    
% trltype: 1=all trials; 2=mismatched trials; 3=matched trials
if trltype==1
    seg_sum = [1,2,3,4,5,6,7,8,9];
elseif trltype==2
    seg_sum = [2,3,4,6,7,8];
elseif trltype==3
    seg_sum = [1,5,9];
end

ind_group = cell(length(seg_sum),2);
for inum = 1:length(seg_sum)
    idx = (label_good==1) & ismember(label_seg,seg_sum(inum));
    ind_group{inum,1} = idx;
    ind_group{inum,2} = sum(idx);
end

n_comb = floor(min(cell2mat(ind_group(:,2)))/kfold);
cvind_group = cell(kfold,length(seg_sum));
for inum = 1:length(seg_sum)
    wholeind_group = find(ind_group{inum,1}==1);
    newind_group = randperm(length(wholeind_group));
    for ifold = 1:kfold
        i_s = (ifold-1)*n_comb + 1;
        i_e = ifold*n_comb; 
        cvind_group{ifold,inum} = wholeind_group(newind_group(i_s:i_e));
    end
end

end