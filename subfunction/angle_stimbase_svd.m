function [cos_pcangle,rev] = angle_stimbase_svd(alldata_run1,alldata_run2,timerange,num_pc)
    
    % averaged across the time
    alldata_run1 = alldata_run1(:,:,timerange);
    alldata_run2 = alldata_run2(:,:,timerange);
    if length(timerange)>1
        alldata_run1 = mean(alldata_run1,3);
        alldata_run2 = mean(alldata_run2,3);
    end
    
    % the first principal angle
    tmp1 = alldata_run1;  % 9*sensor
    tmp1 = tmp1 - mean(tmp1,1);
    [U_1,~,~] = svd((tmp1'*tmp1)/(size(tmp1,1)-1));
    tmp2 = alldata_run2;
    tmp2 = tmp2 - mean(tmp2,1);
    [U_2,~,~] = svd((tmp2'*tmp2)/(size(tmp2,1)-1));
    [~,S,~] = svd(U_1(:,1:num_pc)'*U_2(:,1:num_pc));
    cos_pcangle = S(eye(num_pc)==1);
    clear alldata_run1 alldata_run2

    % ratio of explained variance
    % Li 2023 Curr.Biol.
    rev1 = (std(reshape(tmp1*U_2(:,1:num_pc),[],1))/std(reshape(tmp1*U_1(:,1:num_pc),[],1)))^2;
    rev2 = (std(reshape(tmp2*U_1(:,1:num_pc),[],1))/std(reshape(tmp2*U_2(:,1:num_pc),[],1)))^2;
    rev = mean([rev1,rev2]);

end