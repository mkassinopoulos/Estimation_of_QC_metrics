
clear, clc

NC = 285;  load([ '../../../../Atlases/Gordon_atlas_hcp/Gordon333.mat']); nVect = length(FC_prior_vector); NC_0 = 12;
save_path = ['Export/2019_08_01/Gordon_333/FCC_per_RSN/'];  if exist(save_path,'file')~=7, mkdir(save_path); end

task_list = {'Rest1_LR','Rest1_RL','Rest2_LR','Rest2_RL'};
load('Export/subject_list_390_in_10.mat')

save_path_filename = [save_path,'FCC_GS_WM_200.mat']; 
 
flag_phys = 0;

%%   Optimize pipeline  !!

indCoupl = find(FC_prior_vector==1);
FCC_all = zeros(nScans,1);
FCC_RSN = zeros(nScans,NC_0);

tic
parfor c = 1 : 4*nSubj
    s = ceil(c/4);        run = c - (s-1)*4;
    subject = char(subject_list(s,:));         task = char(task_list(run));
    fprintf('Subject: %s     (%d/%d);   Run: %d/%d    \n',subject,s,nSubj,run,4)
    
    [data, GS, WMpca, CSFpca, FD,movRegr, ~, ~, ~,~,PCAexpl] = load_scan(subject,task,0);
    nComp = size(WMpca,2);           NV = length(GS);             FDmean(c) = mean(FD);
                      
    regr = [ones(NV,1) ,GS, WMpca(:,1:200) ];
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
    
    poolNS = FC_vector; poolNS(indCoupl)=[];
    poolS = FC_vector(indCoupl);
    [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
    FCC_all(c) = a.zval;
    
    FCcontrast_perNetwork = zeros(NC_0,1);
    for r=1:NC_0
        indS = find(FC_prior_perNetwork_vector==r);           poolS = FC_vector(indS);
        [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
        FCcontrast_perNetwork(r) = a.zval;
    end
    FCC_RSN(c,:) = FCcontrast_perNetwork ;
        
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60),  

%  Save stats_metrics ------------------------------

save(save_path_filename,'FCC_all','FCC_RSN')
fprintf('The End \n'),
load chirp,  sound(y,Fs)























    