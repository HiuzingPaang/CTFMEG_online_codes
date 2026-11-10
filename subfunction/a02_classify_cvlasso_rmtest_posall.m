function fold_hit = a02_classify_cvlasso_rmtest_posall(whe_shuffle,setlambda,n_comb,kfold,xall_group,nsensor,gltype)

% gltype: 1=train_global_test_local; 2=train_local_test_global

label_seg = [1,2,3,4,5,6,7,8,9];
label_global = [1,1,1,2,2,2,3,3,3];
label_local = [1,2,3,1,2,3,1,2,3];

if gltype==1
    ytrain_label = label_global;
    ytest_label = label_local;
elseif gltype==2
    ytrain_label = label_local;
    ytest_label = label_global;
end

num = size(xall_group,2);  % xall_group = nan(n_comb*kfold,num,nsensor);
num_classifier = 3;
fold_hit = nan(kfold,1);  % overall
for ifold = 1:kfold
    % define train and test dataset
    i_s = (ifold-1)*n_comb + 1;
    i_e = ifold*n_comb;
    ind_test = false(kfold*n_comb,1);
    ind_test(i_s:i_e,1) = true;
    xtest = reshape(xall_group(ind_test,:,:),n_comb*num,nsensor);
    
    ytest = reshape(repmat(ytest_label,n_comb,1),n_comb*num,1);
    tmp = reshape(repmat(label_seg,n_comb,1),n_comb*num,1);
    xtest = xtest(~ismember(tmp,[1,5,9]),:);
    ytest = ytest(~ismember(tmp,[1,5,9]),1);
    if whe_shuffle==1
        rand_idx = randperm(length(ytest));
        ytest = ytest(rand_idx);
    end
    
    xtrain = reshape(xall_group(~ind_test,:,:),(kfold-1)*n_comb*num,nsensor);
    ytrain = reshape(repmat(ytrain_label,(kfold-1)*n_comb,1),(kfold-1)*n_comb*num,1);
    if whe_shuffle==1
        ytrain = ytrain(randperm(length(ytrain)));
    end

    % lasso, training
    coef_orig = nan(1+nsensor,num_classifier);
    pred_all = nan(n_comb*num - sum(ismember(tmp,[1,5,9])),num_classifier);
    pred = nan(n_comb*num - sum(ismember(tmp,[1,5,9])),1);
    for inum = 1:num_classifier
        % training
        [B,FitInfo] = lassoglm(xtrain,ytrain==inum,'binomial','Alpha',1,'Lambda',setlambda,'Standardize',1);
        B0 = FitInfo.Intercept;
        coef_orig(:,inum) = [B0;B];
        % testing
        pred_all(:,inum) = glmval(coef_orig(:,inum),xtest,'logit');
    end

    for idata = 1:n_comb*num - sum(ismember(tmp,[1,5,9]))
        maxprobind = find(pred_all(idata,:)==max(pred_all(idata,:)));
        if isscalar(maxprobind)
            pred(idata,1) = maxprobind;
        elseif length(maxprobind) > 1
            pred(idata,1) = -1;
        end
    end
    
    fold_hit(ifold,1) = sum(ytest==pred)/length(ytest);
end
clear xall_group xtest xtrain
    
fold_hit = squeeze(mean(fold_hit,1));

end