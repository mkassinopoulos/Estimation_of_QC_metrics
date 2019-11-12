

clear, clc

NC = 285;  load([ 'Atlases/Gordon/Gordon333.mat']); nVect = length(FC_prior_vector);
save_path = ['Export/2019_11_11/Gordon_333/S5_Low_Pass_Filtering/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load('Export/subject_list_390_in_10.mat')

path_name = 'Low_pass_filtering_GS_WM_200' ;    save_path_filename = [save_path,path_name,'.mat'];

%% ====================================

low_pass_list = [ inf, 0.6,0.5, 0.4, 0.3, 0.2, 0.1, 0.08,  0.05, 0.01];
FIX_denoised_flag = 0;
nPipel = length(low_pass_list)*2;
FC_opt_all_vector = zeros(nVect,nScans,nPipel);
FDmean = zeros(nScans,1);
FDDVARS = zeros(nScans,nPipel);
nPipel = size(FC_opt_all_vector,3);
FCC = zeros(nScans,nPipel);
indCoupl = find(FC_prior_vector==1);

TR = 0.72;

tic
parfor c = 1 :  4*nSubj
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [data, GS, WMpca, ~, FD,movRegr] = load_scan(subject,task,0);
    nComp = size(WMpca,2);              NV = length(GS);             FDmean(c) = mean(FD);
    
    regr = [ones(NV,1), GS, WMpca(:,1:200)];
    FC_opt_all_vector_tmp = zeros(nVect,nPipel);
    FDDVARS_tmp = zeros(nPipel,1);
    % LPF before Regression
    
    for model_c = 1:nPipel/2
        f_LPF = low_pass_list(model_c);
        
        data_filt = data; regr_filt = regr;
        if  f_LPF ~= inf
            [filt_b,filt_a] = butter(2,f_LPF*2*TR,'low');
            data_filt = filtfilt(filt_b, filt_a, data);
            regr_filt = filtfilt(filt_b, filt_a, regr);
        end
        
        ROI_data_clean = zeros(size(data));
        for i = 1:NC
            voxel = data_filt(:,i);                B=regr_filt\voxel;   yPred=regr_filt*B;
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
        %         FC_opt_all_vector(:,c,model_c) = FC_vector ;
        
        FC_opt_all_vector_tmp(:,model_c) = FC_vector;
        
        img_diff_col = zeros(size(ROI_data_clean));
        for vox = 1:NC
            voxel = ROI_data_clean(:,vox);
            tmp = diff(voxel);
            img_diff_col(:,vox) = [0;tmp];
        end
        DVARS = rms(img_diff_col,2); DVARS(1) = DVARS(2);
        %         FDDVARS(c,model_c) = corr(FD,DVARS);
        FDDVARS_tmp(model_c) = corr(FD,DVARS);            
    end
    
    
    % LPF after Regression
    
    ROI_data_clean = zeros(size(data));
    for i = 1:NC
        voxel = data(:,i);                B=regr\voxel;   yPred=regr*B;
        ROI_data_clean(:,i)=voxel-yPred;
    end
    
    for model_c = [(nPipel/2 +1) : nPipel]
        f_LPF = low_pass_list(model_c-nPipel/2);
        
        if  f_LPF ~= inf
            [filt_b,filt_a] = butter(2,f_LPF*2*TR,'low');
            ROI_data_clean = filtfilt(filt_b, filt_a, ROI_data_clean);
        end
        FC_tmp = corr(ROI_data_clean);
        
        FC_vector =zeros(nVect,1); k=0;
        for i=1:NC-1
            for j=i+1:NC
                k = k+1;
                FC_vector(k) =FC_tmp(i,j);
            end
        end
        FC_opt_all_vector_tmp(:,model_c) = FC_vector;
        
        img_diff_col = zeros(size(ROI_data_clean));
        for vox = 1:NC
            voxel = ROI_data_clean(:,vox);
            tmp = diff(voxel);
            img_diff_col(:,vox) = [0;tmp];
        end
        DVARS = rms(img_diff_col,2); DVARS(1) = DVARS(2);
        FDDVARS_tmp(model_c) = corr(FD,DVARS);        
    end
    
    for p = 1:nPipel
        FC_opt_all_vector(:,c,p) = FC_opt_all_vector_tmp(:,p);
        
        poolNS = FC_vector; poolNS(indCoupl)=[];
        poolS = FC_vector(indCoupl);
        [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
        FCC(c,p) = a.zval;        
    end        
    FDDVARS(c,:) = FDDVARS_tmp;
    
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60)


%%  Estimate and Save QC metrics      ------------

Estimate_QC_metrics_for_all_groups
save(save_path_filename,'stats_metrics_all_groups')
load chirp,  sound(y,Fs),














