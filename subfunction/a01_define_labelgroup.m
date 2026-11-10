function [num_group,label_group] = a01_define_labelgroup(classtype,numrun,behmeglog,stimcellpath,subid)

if strcmp(classtype,'global')
    num_group = 3;
    label_group = repmat([1,1,1,2,2,2,3,3,3],1,9*numrun)';
elseif strcmp(classtype,'local')
    num_group = 3;
    label_group = repmat([1,2,3,1,2,3,1,2,3],1,9*numrun)';
elseif strcmp(classtype,'seg')
    num_group = 9;
    label_group = repmat(1:9,1,9*numrun)';
elseif strcmp(classtype,'scene')
    num_group = 3;
    % get the label_group of scene
    load([stimcellpath,'/sub',num2str(subid),'WholeListStimcell.mat'],'stimcell');
    vidind = find(strcmp(stimcell(1,:),{'LeaSce1'}));
    scenelabel_pre = cell2mat(stimcell(2:end,vidind:vidind+2));
    scenelabel_post = [repmat(scenelabel_pre(:,1),1,3),repmat(scenelabel_pre(:,2),1,3),repmat(scenelabel_pre(:,3),1,3)];
    scenelabel = nan(numrun*9,9);
    for irun = 1:numrun
        is0 = (irun-1)*9+1;
        ie0 = irun*9;
        is1 = (behmeglog(1,irun)-1)*9+1;
        ie1 = behmeglog(1,irun)*9;
        scenelabel(is0:ie0,:) = scenelabel_post(is1:ie1,:);
    end
    label_group = reshape(scenelabel',[],1);
elseif strcmp(classtype,'animal')
    num_group = 15;
    % get the label_group of animal
    load([stimcellpath,'/sub',num2str(subid),'WholeListStimcell.mat'],'stimcell');
    sceind = find(strcmp(stimcell(1,:),{'LeaSce1'}));
    aniind = find(strcmp(stimcell(1,:),{'LeaAni1'}));
    scelabel_pre = cell2mat(stimcell(2:end,sceind:sceind+2));
    anilabel_pre = cell2mat(stimcell(2:end,aniind:aniind+2));
    animal_pre = (scelabel_pre-1)*5+anilabel_pre;
    animal_post = [repmat(animal_pre(:,1),1,3),repmat(animal_pre(:,2),1,3),repmat(animal_pre(:,3),1,3)];
    animallabel = nan(numrun*9,9);
    for irun = 1:numrun
        is0 = (irun-1)*9+1;
        ie0 = irun*9;
        is1 = (behmeglog(1,irun)-1)*9+1;
        ie1 = behmeglog(1,irun)*9;
        animallabel(is0:ie0,:) = animal_post(is1:ie1,:);
    end
    label_group = reshape(animallabel',[],1);
end

end