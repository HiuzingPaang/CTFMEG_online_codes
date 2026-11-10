function fold_hit = a01_classify_cvlasso_pos23(whe_shuffle,setlambda,n_comb,kfold,xall_group,nsensor,trltype,gltype)

% gltype: 1=global; 2=local
% trltype: 1=all trials; 2=mismatched trials

if trltype==1
    label_seg = [1,2,3,4,5,6,7,8,9];
    label_global = [1,1,1,2,2,2,3,3,3];
    label_local = [1,2,3,1,2,3,1,2,3];
elseif trltype==2
    label_seg = [2,3,4,6,7,8];
    label_global = [1,1,2,2,3,3];
    label_local = [2,3,1,3,1,2];
end

if gltype==1
    y_label = label_global;
elseif gltype==2
    y_label = label_local;
end

y_label = y_label(y_label~=1);
label_seg = label_seg(y_label~=1);
xall_group = xall_group(:,y_label~=1,:);
num = sum(y_label~=1);
fold_hit = nan(kfold,3+num);  % overall, positions 2-3, serial positions
for ifold = 1:kfold
    % define train and test dataset
    i_s = (ifold-1)*n_comb + 1;
    i_e = ifold*n_comb;
    ind_test = false(kfold*n_comb,1);
    ind_test(i_s:i_e,1) = true;
    xtest = reshape(xall_group(ind_test,:,:),n_comb*num,nsensor);
    
    ytest = reshape(repmat(y_label,n_comb,1),n_comb*num,1);
    ytest_seg_label = reshape(repmat(label_seg,n_comb,1),n_comb*num,1);
    if whe_shuffle==1
        rand_idx = randperm(length(ytest));
        ytest = ytest(rand_idx);
        ytest_seg_label = ytest_seg_label(rand_idx);
    end
    
    xtrain = reshape(xall_group(~ind_test,:,:),(kfold-1)*n_comb*num,nsensor);
    ytrain = reshape(repmat(y_label,(kfold-1)*n_comb,1),(kfold-1)*n_comb*num,1);
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
    for inum = 1:num
        fold_hit(ifold,3+inum) = sum((ytest_seg_label==label_seg(inum)) & (ytest==pred))/sum(ytest_seg_label==label_seg(inum));
    end
    
end
clear xall_group xtest xtrain
    
fold_hit = squeeze(mean(fold_hit,1));

end