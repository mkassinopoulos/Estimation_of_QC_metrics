
clear, clc

NC = 285;  load([ 'Atlases/Gordon/Gordon333.mat']); nVect = length(FC_prior_vector);
save_path = ['Export/2019_11_11/Gordon_333/S4_Scrubbing/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load('Export/subject_list_390_in_10.mat')

flag_FD_1_DVARS_2 = 1;    %   Set it to 1 or 2 to consider FD or DVARS as a framewise index of quality to flag motion-contaminated volumes
path_name = 'GS_WM_200_FD' ;    save_path_filename = [save_path,path_name,'.mat'];

%% ====================================

if  flag_FD_1_DVARS_2 == 1
    FD_thresh_list = [0.1:0.05:1.0,inf]   % FD
else
    FD_thresh_list = [0.5:0.5:20,inf]        % DVARS
end

nPipel = length(FD_thresh_list);
FC_opt_all_vector = zeros(nVect,nScans,nPipel);
FDmean = zeros(nScans,1);
FDDVARS = zeros(nScans,nPipel);
nPipel = size(FC_opt_all_vector,3);

tic
parfor c = 1 : 4*nSubj
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [data, GS, ~, WMpca, ~, ~,  FD,movRegr, ~, ~,~,DVARS_WB] = load_scan(subject,task,0);
    nComp = size(WMpca,2);              NV = length(GS);             FDmean(c) = mean(FD);
    
    for model_c = 1:nPipel   % par
        
        if flag_FD_1_DVARS_2 == 1
            ind =find(FD<FD_thresh_list(model_c));    % only FD
        else
            ind_rm =   find(isoutlier(DVARS_WB,'median','ThresholdFactor', FD_thresh_list(model_c)));  ind=1:NV; ind(ind_rm) = [];    % only DVARS
        end
        
        regr = [ones(NV,1),GS, WMpca(:,1:200)];
        regr = regr(ind,:);
        data_tmp = data(ind,:);
        ROI_data_clean = zeros(length(ind),NC);
        for i = 1:NC
            voxel = data_tmp(:,i);                B=regr\voxel;   yPred=regr*B;
            ROI_data_clean(:,i)=voxel-yPred;
        end
        FC_tmp = corr(ROI_data_clean);
        
        FC_vector =zeros(nVect,1); k=0;
        for i=1:NC-1
            for j=i+1:NC
                k = k+1;
                FC_vector(k) =FC_tmp(i,j);
            end
        end
        FC_opt_all_vector(:,c,model_c) = FC_vector ;
        
        img_diff_col = zeros(size(ROI_data_clean));
        for vox = 1:NC
            voxel = ROI_data_clean(:,vox);
            tmp = diff(voxel);
            img_diff_col(:,vox) = [0;tmp];
        end
        DVARS = rms(img_diff_col,2); DVARS(1) = DVARS(2);
        FDDVARS(c,model_c) = corr(FD(ind),DVARS);
    end
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60)


%%  Estimate and Save QC metrics      ------------

Estimate_QC_metrics_for_all_groups
save(save_path_filename,'stats_metrics_all_groups')
load chirp,  sound(y,Fs),
















