

clear, clc

NC = 285;  load([ 'Atlases/Gordon/Gordon333.mat']); nVect = length(FC_prior_vector);
save_path = ['Export/2019_11_11/Gordon_333/S3_PCA_vs_ModelBased_NCTs/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load('Export/subject_list_390_in_10.mat')

path_name = 'PNC_models_64' ;         save_path_filename = [save_path,'PCA_',path_name,'.mat'];

%% ====================================


load('Struct_64.mat'),   nPipel = size(struct_64,1);
FC_opt_all_vector = zeros(nVect,nScans,nPipel);
FDmean = zeros(nScans,1);
FDDVARS = zeros(nScans,nPipel);

TR = 0.72;

tic
parfor c = 1 : 4*nSubj
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [data, GS,~, WMpca, ~, ~,  FD,movRegr, RETR_RespRegr, RETR_CardRegr, SLFOs] = load_scan(subject,task,0);
    nComp = size(WMpca,2);              NV = length(GS);             FDmean(c) = mean(FD);
    
    TBpca = WMpca(:,1:200);
    
    for model_c = 1:nPipel
        DM = struct_64(model_c,:);
        regr = [ones(NV,1)];
        if DM(1)==1,   regr = [regr, RETR_CardRegr]; end
        if DM(2)==1,   regr = [regr, RETR_RespRegr]; end
        if DM(3)==1,   regr = [regr, movRegr]; end
        if DM(4)==1,   regr = [regr, SLFOs]; end
        if DM(5)==1,   regr = [regr, GS]; end
        if DM(6)==1,   regr = [regr, TBpca]; end
        
        ROI_data_clean = zeros(size(data));
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
        FC_opt_all_vector(:,c,model_c) = FC_vector ;
        
        img_diff_col = zeros(size(ROI_data_clean));
        for vox = 1:NC
            voxel = ROI_data_clean(:,vox);
            tmp = diff(voxel);
            img_diff_col(:,vox) = [0;tmp];
        end
        DVARS = rms(img_diff_col,2); DVARS(1) = DVARS(2);
        FDDVARS(c,model_c) = corr(FD,DVARS);
    end
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60), load chirp,  sound(y,Fs)


%%  Estimate and Save QC metrics      ------------

Estimate_QC_metrics_for_all_groups
save(save_path_filename,'stats_metrics_all_groups')
load chirp,  sound(y,Fs),










