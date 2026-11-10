function fold_hit = a02_classify_cvlasso_pos23(whe_shuffle,setlambda,n_comb,kfold,xall_group,nsensor,gltype,trltype)

% gltype: 1=train_global_test_local; 2=train_local_test_global
% trltype: 1=all trials; 2=mismatched trials; 3=matched trials

if trltype==1
    label_seg = [1,2,3,4,5,6,7,8,9];
    label_global = [1,1,1,2,2,2,3,3,3];
    label_local = [1,2,3,1,2,3,1,2,3];
elseif trltype==2
    label_seg = [2,3,4,6,7,8];
    label_global = [1,1,2,2,3,3];
    label_local = [2,3,1,3,1,2];
elseif trltype==3
    label_seg = [1,5,9];
    label_global = [1,2,3];
    label_local = [1,2,3];
end

if gltype==1
    ytrain_label = label_global;
    ytest_label = label_local;
elseif gltype==2
    ytrain_label = label_local;
    ytest_label = label_global;
end
num_seg_pos = label_seg(ytest_label~=1);

num = size(xall_group,2);  % xall_group = nan(n_comb*kfold,num,nsensor);
fold_hit = nan(kfold,3+length(num_seg_pos));  % overall, positions 2-3, serial positions
for ifold = 1:kfold
    % define train and test dataset
    i_s = (ifold-1)*n_comb + 1;
    i_e = ifold*n_comb;
    ind_test = false(kfold*n_comb,1);
    ind_test(i_s:i_e,1) = true;
    xtest = reshape(xall_group(ind_test,:,:),n_comb*num,nsensor);
    
    ytest = reshape(repmat(ytest_label,n_comb,1),n_comb*num,1);
    ytest_seg_label = reshape(repmat(label_seg,n_comb,1),n_comb*num,1);
    xtest = xtest(ismember(ytest,[2,3]),:);
    ytest_seg_label = ytest_seg_label(ismember(ytest,[2,3]),1);
    ytest = ytest(ismember(ytest,[2,3]),1);
    if whe_shuffle==1
        rand_idx = randperm(length(ytest));
        ytest = ytest(rand_idx);
        ytest_seg_label = ytest_seg_label(rand_idx);
    end

    xtrain = reshape(xall_group(~ind_test,:,:),(kfold-1)*n_comb*num,nsensor);
    ytrain = reshape(repmat(ytrain_label,(kfold-1)*n_comb,1),(kfold-1)*n_comb*num,1);
    xtrain = xtrain(ismember(ytrain,[2,3]),:);
    ytrain = ytrain(ismember(ytrain,[2,3]),1);
    if whe_shuffle==1
        ytrain = ytrain(randperm(length(ytrain)));
    end
    
    % lasso, training
    [B,FitInfo] = lassoglm(xtrain,ytrain==3,'binomial','Alpha',1,'Lambda',setlambda,'Standardize',1);
    B0 = FitInfo.Intercept;
    % testing
    pred = double(glmval([B0;B],xtest,'logit') > 0.5) + 2;
        
    fold_hit(ifold,1) = sum(ytest==pred)/length(ytest);
    for inum = 2:3
        fold_hit(ifold,inum) = sum((ytest==inum) & (ytest==pred))/sum(ytest==inum);
    end
    for inum = 1:length(num_seg_pos)
        fold_hit(ifold,3+inum) = sum((ytest_seg_label==num_seg_pos(inum)) & (ytest==pred))/sum(ytest_seg_label==num_seg_pos(inum));
    end
end
clear xall_group xtest xtrain
    
fold_hit = squeeze(mean(fold_hit,1));

end