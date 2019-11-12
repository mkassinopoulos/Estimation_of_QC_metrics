
function stats  = Estimate_QC_metrics(nPipel_f, nScans_f, nSubj_f,FC_opt_all_vector_f,FDmean_f)


%%  Get first scan from each subject

ind = 1:4:nScans_f;
FC_all_vector_subj = FC_opt_all_vector_f(:,ind,:);
FDmean_subj = FDmean_f(ind);

%%  Estimate: FCC, FCC_per_network, FD_FCC  ----------

NC = 285;  load([ 'Atlases/Gordon/Gordon333.mat']); nVect = length(FC_prior_vector); NC_0 = 12;

FC_prior_perNetwork_vector = FC_prior_perNetwork_vector;
indCoupl = find(FC_prior_vector==1);
FCC_all = zeros(nSubj_f,nPipel_f,NC_0+1);

parfor c = 1 : nSubj_f
    FCC_tmp = zeros(nPipel_f,NC_0+1);
    for p = 1 : nPipel_f
        FC_vector = squeeze(FC_all_vector_subj(:,c,p));
        poolNS = FC_vector; poolNS(indCoupl)=[];
        poolS = FC_vector(indCoupl);
        [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
        FCC_tmp(p,1) = a.zval;
        
        FCcontrast_perNetwork = zeros(NC_0,1);
        for r=1:NC_0
            indS = find(FC_prior_perNetwork_vector==r);           poolS = FC_vector(indS);
            [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
            FCcontrast_perNetwork(r) = a.zval;
        end
        FCC_tmp(p,2:NC_0+1) = FCcontrast_perNetwork;
    end
    FCC_all(c,:,:) = FCC_tmp ;
end

tmp = zeros(NC_0,nPipel_f);
for i = 1:NC_0
    tmp(i,:) = mean(squeeze(FCC_all(:,:,i+1)),1);
end
stats.FCCperN = tmp';

FCC_subj = squeeze(FCC_all(:,:,1));
stats.FCCall = mean(FCC_subj,1)';

tmp = zeros(nPipel_f,1);
for p = 1:nPipel_f
    x = FCC_subj(:,p);
    tmp(p) = corr(x,FDmean_subj);
end
stats.FD_FCC = tmp;

%% Estimate:   FDFC_median, FCFD_distance, FD-MFC

r_FDFCmedian = zeros(nPipel_f,1);
r_FD_MFC = zeros(nPipel_f,1);
FDFC_distance = zeros(nPipel_f,1);

for p = 1:nPipel_f
    FC_vector = squeeze(FC_all_vector_subj(:,:,p));
    FDFC_vector = zeros(nVect,1);
    parfor i = 1:nVect
        FC_vector_edge = FC_vector(i,:)';
        FDFC_vector(i) = corr(FC_vector_edge, FDmean_subj);
    end
    r_FDFCmedian(p) = median(abs(FDFC_vector(:)));
    FCmean = mean(FC_vector,1)';
    r_FD_MFC(p) = corr(FDmean_subj, FCmean);
    FDFC_distance(p) = corr(dist_vect,FDFC_vector);
end

stats.FDFCmedian = r_FDFCmedian;
stats.FDFCmean_v2 = r_FD_MFC;
stats.FDFC_distance = FDFC_distance;


%%   Intraclass correlation for within-subject reliability
%  Estimate: Median ICC (MICCC),  ICC contrast (ICCC)

M_ICC = zeros(nVect,nPipel_f,nSubj_f,4);
for subject=1:nSubj_f
    for sc=1:4
        for xInd = 1:nPipel_f
            for i = 1:nVect
                M_ICC(i,xInd,subject,sc) = FC_opt_all_vector_f(i,(subject-1)*4+sc,xInd);
            end
        end
    end
end

ICCstats = zeros(nVect,nPipel_f);
parfor i = 1:nVect
    for xInd = 1:nPipel_f
        M = squeeze(M_ICC(i,xInd,:,:));
        ICCstats(i,xInd) = ICC(M,'1-k',0.05);
    end
end
stats_median = median(ICCstats);

ICCcontrast = zeros(nPipel_f,1);
for k = 1:nPipel_f
    FC_vector = ICCstats(:,k);
    poolNS = FC_vector; poolNS(indCoupl)=[];
    poolS = FC_vector(indCoupl);
    [ttest_p,ttest_h,a] = ranksum(poolS,poolNS,'Tail','right');
    ICCcontrast(k) = a.zval  ;
end

stats.ICCmedian =  stats_median';
stats.ICCcontrast =  ICCcontrast;


%%  ------------------------












