function [fold_coef_orig,fold_coef_corr,fold_hit,pred_trlacc,pred_trlacc_all] = a01_classify_cvlasso_posall(whe_shuffle,setlambda,n_comb,kfold,xall_group,nsensor,xall_group_trlind,trltype,gltype)

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

num_classifier = 3;
num = length(label_seg);
% create the aimed matrix
fold_coef_orig = nan(kfold,1+nsensor,num_classifier);
fold_coef_corr = nan(kfold,nsensor,num_classifier);
fold_hit = nan(kfold,1+3+num);  % overall, positions 1-3, serial positions
pred_trlacc = nan(kfold*n_comb*num,4);  % 1=trial index, 2=trial acc, 3=true label, 4=predicted label(-1,1,2,3,...)
pred_trlacc_all = nan(kfold*n_comb*num,2+num_classifier);  % 1=trial index, 2=true label, 3=predicted probability from all classifiers

for ifold = 1:kfold
    % define train and test dataset
    i_s = (ifold-1)*n_comb + 1;
    i_e = ifold*n_comb;
    ind_test = false(kfold*n_comb,1);
    ind_test(i_s:i_e,1) = true;
    xtest = reshape(xall_group(ind_test,:,:),n_comb*num,nsensor);
    
    test_trlind = reshape(xall_group_trlind(ind_test,:),n_comb*num,1);    
    ii_s = (ifold-1)*n_comb*num + 1;
    ii_e = ifold*n_comb*num;
    pred_trlacc(ii_s:ii_e,1) = test_trlind;
    pred_trlacc_all(ii_s:ii_e,1) = test_trlind;
    
    ytest = reshape(repmat(y_label,n_comb,1),n_comb*num,1);
    ytest_seg_label = reshape(repmat(label_seg,n_comb,1),n_comb*num,1);
    if whe_shuffle==1
        rand_idx = randperm(length(ytest));
        ytest = ytest(rand_idx);
        ytest_seg_label = ytest_seg_label(rand_idx);
    end
    
    pred_trlacc(ii_s:ii_e,3) = ytest;
    pred_trlacc_all(ii_s:ii_e,2) = ytest;
    
    xtrain = reshape(xall_group(~ind_test,:,:),(kfold-1)*n_comb*num,nsensor);
    ytrain = reshape(repmat(y_label,(kfold-1)*n_comb,1),(kfold-1)*n_comb*num,1);
    if whe_shuffle==1
        ytrain = ytrain(randperm(length(ytrain)));
    end

    % lasso
    coef_orig = nan(1+nsensor,num_classifier);
    coef_corr = nan(nsensor,num_classifier);
    pred_all = nan(n_comb*num,num_classifier);
    pred = nan(n_comb*num,1);
    for inum = 1:num_classifier
        % training
        [B,FitInfo] = lassoglm(xtrain,ytrain==inum,'binomial','Alpha',1,'Lambda',setlambda,'Standardize',1);
        B0 = FitInfo.Intercept;
        coef_orig(:,inum) = [B0;B];
        % testing
        pred_all(:,inum) = glmval(coef_orig(:,inum),xtest,'logit');
        % reconstruct the activation pattern
        coef_corr(:,inum) = cov(xtrain)*coef_orig(2:end,inum);
    end

    for idata = 1:n_comb*num
        maxprobind = find(pred_all(idata,:)==max(pred_all(idata,:)));
        if isscalar(maxprobind)
            pred(idata,1) = maxprobind;
        elseif length(maxprobind) > 1
            pred(idata,1) = -1;
        end
    end
    
    fold_hit(ifold,1) = sum(ytest==pred)/length(ytest);
    for inum = 1:num_classifier
        fold_hit(ifold,1+inum) = sum((ytest==inum) & (ytest==pred))/sum(ytest==inum);
    end
    for inum = 1:num
        fold_hit(ifold,4+inum) = sum((ytest_seg_label==label_seg(inum)) & (ytest==pred))/sum(ytest_seg_label==label_seg(inum));
    end

    pred_trlacc(ii_s:ii_e,2) = double(ytest==pred);
    pred_trlacc(ii_s:ii_e,4) = pred;
    pred_trlacc_all(ii_s:ii_e,3:end) = pred_all;
    fold_coef_orig(ifold,:,:) = coef_orig;
    fold_coef_corr(ifold,:,:) = coef_corr;
end
clear xall_group xtest xtrain
    
fold_coef_orig = squeeze(mean(fold_coef_orig,1));
fold_coef_corr = squeeze(mean(fold_coef_corr,1));
fold_hit = squeeze(mean(fold_hit,1));

end