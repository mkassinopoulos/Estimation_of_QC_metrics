%%  Estimate QC metrics

addpath('Libraries/ICC')

nSubj_tmp = size(subjects_in_10,1);
N_SS = size(subjects_in_10,2);
nScans_tmp = nSubj_tmp*4;

% create structure for stats_metrics with random values
stats_metrics =  Estimate_QC_metrics(1,  nScans_tmp, nSubj_tmp, FC_opt_all_vector(:,1:nScans_tmp,1), FDmean(1:nScans_tmp));
stats_metrics.FDDVARS =  0;     stats_metrics.FD_FDDVARS =  0;
clear stats_metrics_all_groups
for i = 1:N_SS,         stats_metrics_all_groups(i) = stats_metrics;    end

parfor i = 1:N_SS      %   in case you have issues with RAM replace parfor with for  (i.e., don't run the loop in parallel)
    fprintf('Sample: %d/%d   \n', i, N_SS),
    ind = subjects_in_10(:,i);
    ind_scans = [1:4]'+(ind'-1)*4; ind_scans = ind_scans(:);
    ind_scans_subj = ind_scans(1:4:end);
    
    FDmean_tmp = FDmean(ind_scans);
    FDmean_subj = FDmean(ind_scans_subj);
    FDDVARS_subj = FDDVARS(ind_scans_subj,:);
    
    FC_opt_all_vector_tmp = FC_opt_all_vector(:,ind_scans,:);
    nScans_tmp = length(ind_scans); nSubj_tmp = nScans_tmp/4;
    
    stats_metrics =  Estimate_QC_metrics(nPipel,  nScans_tmp, nSubj_tmp, FC_opt_all_vector_tmp, FDmean_tmp);
    stats_metrics.FDDVARS =  mean(FDDVARS_subj)';
    
    FD_FDDVARS = zeros(nPipel,1);
    for p = 1:nPipel
        FD_FDDVARS(p) = corr(FDmean_subj, FDDVARS_subj(:,p));
    end
    stats_metrics.FD_FDDVARS =  FD_FDDVARS;
    stats_metrics_all_groups(i) = stats_metrics;
end
fprintf('Time elapsed (minutes): %3.1f  \n', toc/60),