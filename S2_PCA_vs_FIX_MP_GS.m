
clear, clc

NC = 285;  load( 'Atlases/Gordon/Gordon333.mat'); nVect = length(FC_prior_vector);
save_path = ['Export/2019_11_11/Gordon_333/S2_only_one_pipeline/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load('Export/subject_list_390_in_10.mat')

path_name = 'GS_WM_200'

save_path_filename = [save_path,'PCA_',path_name,'.mat']; save_path_filename_FC = [save_path,'FC_',path_name,'.mat'];

%%   Optimize pipeline  !!

FC_opt_all_vector = zeros(nVect,nScans);
FDmean = zeros(nScans,1);
FDDVARS = zeros(nScans,1);
FCC = zeros(nScans,1);
indCoupl = find(FC_prior_vector==1);

tic
parfor c = 1 : nScans
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [data, GS, ~, WMpca, ~, CSFpca, FD] = load_scan(subject,task,0);
    nComp = size(WMpca,2);           NV = length(GS);             FDmean(c) = mean(FD);
    
    regr = [ones(NV,1) , GS, WMpca(:,1:200) ];
    ROI_data_clean = zeros(NV,NC);
    for i = 1:NC
        voxel = data(:,i);                B=regr\voxel;   yPred=regr*B;
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
    FC_opt_all_vector(:,c) = FC_vector ;
    
    img_diff_col = zeros(size(ROI_data_clean));
    for vox = 1:NC
        voxel = ROI_data_clean(:,vox);
        tmp = diff(voxel);
        img_diff_col(:,vox) = [0;tmp];
    end
    DVARS = rms(img_diff_col,2); DVARS(1) = DVARS(2);
    FDDVARS(c) = corr(FD,DVARS);
    
    poolNS = FC_vector; poolNS(indCoupl)=[];    poolS = FC_vector(indCoupl);
    [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
    FCC(c) = a.zval;
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60),  load chirp,  sound(y,Fs)

save(save_path_filename_FC,'FC_all_vector','FDDVARS','FCC','FDmean')   % add '-v7.3' if it gives errors

%%  Estimate and Save QC metrics      ------------

Estimate_QC_metrics_for_all_groups
save(save_path_filename,'stats_metrics_all_groups')
load chirp,  sound(y,Fs),










