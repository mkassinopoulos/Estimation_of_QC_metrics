

clear, clc

NC = 285;  load([ 'Atlases/Gordon/Gordon333.mat']); nVect = length(FC_prior_vector);
save_path = ['Export/2019_11_11/Gordon_333/S2_pipelines_from_papers/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load('Export/subject_list_390_in_10.mat')

path_name = 'PCA_pipelines_from_papers' ;         save_path_filename = [save_path,path_name,'.mat'];


%% ====================================

nPipel = 15;
FC_opt_all_vector = zeros(nVect,nScans,nPipel);
FDmean = zeros(nScans,1);
FDDVARS = zeros(nScans,nPipel);
indCoupl = find(FC_prior_vector==1);
FCC = zeros(nScans,nPipel);

tic
parfor c = 1 : 4*nSubj
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [data, GS, WMm, WMpca, CSFm, CSFpca, FD,movRegr,~,~,~,~,PCAexpl] = load_scan(subject,task,0);
    nComp = size(WMpca,2);              NV = length(GS);             FDmean(c) = mean(FD);
    
    GS_der = [0; diff(GS)];     WM_der = [0; diff(WMm)];     CSF_der = [0; diff(CSFm)];
    GS_2 = GS.^2;      WM_2 = WMm.^2;       CSF_2 = CSFm.^2;
    GS_der_2 = GS_der.^2;      WM_der_2 = WM_der.^2;       CSF_der_2 = CSF_der.^2;
    
    WMexpl = cumsum(PCAexpl(:,1));   
    CSFexpl = cumsum(PCAexpl(:,2));  
    
    %%  =========================================
        
    % P1: MP - 6 par
    regr_1 = [ones(NV,1), movRegr(:,1:6)];
    
    % P2: MP - 12 par
    regr_2 = [ones(NV,1), movRegr(:,1:12)];
    
    % P3: MP - 24 par
    regr_3 = [ones(NV,1), movRegr(:,1:24)];
            
    % P5: MP_12 + GS
    regr_5 = [ones(NV,1), GS, movRegr(:,1:12) ];
    
    % P6: WM_5, CSF_5
    regr_6 = [ones(NV,1) , WMpca(:,1:5), CSFpca(:,1:5)];
    
    % P7: 12 MP, WM, CSF
    regr_7 = [ones(NV,1) , WMm, CSFm, movRegr(:,1:12) ];
    
    % P8: 12 MP, WM, CSF, GS
    regr_8 = [ones(NV,1) , WMm, CSFm, movRegr(:,1:12),  GS];
    
    % P9: 30-par:  [GS, WM, CSF]x[TS, derivatives], 24 MP
    regr_9 = [ones(NV,1), movRegr(:,1:24),  GS, WMm, CSFm,  GS_der, WM_der, CSF_der];
    
    % P10: 36-par: [6 MP, GS, WM, CSF] x [TS, derivatives, quadratic and their derivatives
    regr_10 = [ones(NV,1), movRegr(:,1:24),  GS, WMm, CSFm,  GS_der, WM_der, CSF_der, GS_2, WM_2, CSF_2, GS_der_2,  WM_der_2, CSF_der_2];
            
    % P14:   50% of var in WM and CSF
    [~, loc] = min(abs(WMexpl-50)); ind_comp_WM = 1:min(loc,nComp);   % PCAexplained
    [~, loc] = min(abs(CSFexpl-50)); ind_comp_CSF = 1:min(loc,nComp);   % PCAexplained
    regr_14 = [ones(NV,1) , WMpca(:,ind_comp_WM),  CSFpca(:,ind_comp_CSF) ];    % size(regr_16)
    
    % P15:   30% of var in WM
    [~, loc] = min(abs(WMexpl-30)); ind_comp_WM = 1:min(loc,nComp);   % PCAexplained
    regr_15 = [ones(NV,1) , GS,  WMpca(:,ind_comp_WM) ];    % size(regr_16)
    
    % P16:   35% of var in WM
    [~, loc] = min(abs(WMexpl-35)); ind_comp_WM = 1:min(loc,nComp);   % PCAexplained
    regr_16 = [ones(NV,1) , GS,  WMpca(:,ind_comp_WM) ];    % size(regr_16)
    
    % P17:   40% of var in WM
    [~, loc] = min(abs(WMexpl-40)); ind_comp_WM = 1:min(loc,nComp);   % PCAexplained
    regr_17 = [ones(NV,1) , GS, WMpca(:,ind_comp_WM) ];    % size(regr_16)
    
    % P18:   45% of var in WM
    [~, loc] = min(abs(WMexpl-45)); ind_comp_WM = 1:min(loc,nComp);   % PCAexplained
    regr_18 = [ones(NV,1) , GS,  WMpca(:,ind_comp_WM) ];    % size(regr_16)
    
    % P19:   50% of var in WM
    [~, loc] = min(abs(WMexpl-50)); ind_comp_WM = 1:min(loc,nComp);   % PCAexplained
    regr_19 = [ones(NV,1) , GS, WMpca(:,ind_comp_WM) ];    % size(regr_16)
    
    regr_list = [1:3, 5:10, 14:19];
        regr_all = {regr_1, regr_2, regr_3, regr_5, regr_6,  regr_7, regr_8, regr_9, regr_10, ...
                         regr_14, regr_15, regr_16, regr_17, regr_18, regr_19};

    %%  =========================================
    
    for model_c = 1  : nPipel        
        regr = regr_all{model_c};         
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
        
        poolNS = FC_vector; poolNS(indCoupl)=[];
        poolS = FC_vector(indCoupl);
        [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
        FCC(c,model_c) = a.zval;               
    end      
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60), load chirp,  sound(y,Fs)


%%  Estimate and Save QC metrics      ------------

Estimate_QC_metrics_for_all_groups
save(save_path_filename,'stats_metrics_all_groups')
load chirp,  sound(y,Fs),













