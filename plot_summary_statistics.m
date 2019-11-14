

%%  Choose section (group of files)

close all, clear, clc

%  Sections:
%  1:  Plot QC metrics for varying number of WM and CSF denoising
%  2:  Plot QC metrics for pipelines used in the literature (see Table 1 in
%  Kassinopoulos and Mitsis, bioRxiv) for a description of the pipelines  examined)
%  3:  Plot QC metrics for model-based NCTs (motion parameters, cardiac and
%  respiratory regressors (3rd order RETROICOR)
%  4:  Plot QC metrics for the effect of scrubbing (FD and DVARS)
%  5:  Plot QC metrics for the effect of low-pass filtering

% Choose section

Section = 1

% ----------------------------------------------
Dir_G333_S1order = 'Export\2019_11_11\Gordon_333\S1_PCA_order\';
Dir_G333_S2 = 'Export\2019_11_11\Gordon_333\S2_only_one_pipeline\';
Dir_G333_S2b = 'Export\2019_11_11\Gordon_333\S2_pipelines_from_literature\';
Dir_G333_S3 = 'Export\2019_11_11\Gordon_333\S3_PCA_vs_ModelBased_NCTs\';
Dir_G333_S4 = 'Export\2019_11_11\Gordon_333\S4_Censoring\';
Dir_G333_S5 = 'Export\2019_11_11\Gordon_333\S5_Low_Pass_Filtering\';

N_RSN = 12;
stats_metrics_all = [];

%    -----------------------------------------------------

if Section == 1
    k = 1;  stats(k) = load([Dir_G333_S1order,'PCA_WM_order'])  ;
    k = k+1;  stats(k) = load([Dir_G333_S1order,'PCA_CSF_order'])  ;
    k = k+1;  stats(k) = load([Dir_G333_S1order,'PCA_GS_WM_order'])  ;
    k = k+1;  stats(k) = load([Dir_G333_S1order,'PCA_GS_CSF_order'])  ;
end

if Section == 2
    k=1;    stats(k) = load([Dir_G333_S2,'Raw'])  ;
    k = k+1; stats(k) = load([Dir_G333_S2,'GS'])  ;
    k = k+1;   stats(k) = load([Dir_G333_S2,'FIX'])  ;
    k = k+1;  stats(k) = load([Dir_G333_S2,'FIX_GS'])  ;
    k = k+1;  stats(k) = load([Dir_G333_S2,'PCA_GS_WM_200']) ;
    k = k+1;  stats(k) = load([Dir_G333_S2b,'QC_for_pipelines_from_literature.mat'])  ;
end

if Section == 3
    k = 1;  stats(k) = load([Dir_G333_S3,'PCA_vs_ModelBased_NCTs'])  ;
end

if Section == 4   %
    k = 1;  stats(k) = load([Dir_G333_S4,'GS_WM_200_FD'])  ;
    k = k+1;  stats(k) = load([Dir_G333_S4,'GS_WM_200_DVARS'])  ;
end

if Section == 5   %
    k = 1;  stats(k) = load([Dir_G333_S5,'GS_WM_200_LPF_bef_and_after_regr'])  ;
end

%      -----------------------------------------------------------

Dir_x = Dir_G333_S2;

k = k+1;    stats(k) = load([Dir_G333_S2,'Raw'])  ;
k = k+1; stats(k) = load([Dir_G333_S2,'GS'])  ;
k = k+1; stats(k) = load([Dir_G333_S2,'FIX'])  ;
k = k+1; stats(k) = load([Dir_G333_S2,'FIX_GS'])  ;

N_GP = size(stats,2); N_groups = zeros(N_GP,1);
for i = 1:N_GP
    N_groups(i) =  length(stats(i).stats_metrics_all_groups(1).FCCall);
end

clear ranges
for i = 1:length(N_groups)
    if i ==1
        ranges(i,1) = 1; ranges(i,2) = N_groups(i);
    else
        ranges(i,1) = ranges(i-1,1)+N_groups(i-1);
        ranges(i,2) = ranges(i,1) + N_groups(i)-1;
    end
end
nModels = sum(N_groups)-4;
ind_raw = nModels +1;

%  reformat matrices  ---------------------------------

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FCCall;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
FCC = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FD_FCC;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
FD_FCC = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).ICCmedian;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
ICCMed = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).ICCcontrast;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
ICCcontr = x;

%  -------------------------------------------


x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FDDVARS;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
xFDDVARS = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FD_FDDVARS;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
xFD_FDDVARS = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FDFCmedian;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
xFDFCmed = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FDFCmean_v2;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
xFDFCM2 = x;

x = [];
for gp = 1:N_GP
    tmp = [];
    for gs = 1:10
        tmp_1 = stats(gp).stats_metrics_all_groups(gs).FDFC_distance;
        tmp = [tmp; tmp_1(:)'];
    end
    x = [x, tmp];
end
xFDFCdist = x;

%  Extract normalised QC metrics


x = FCC;    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
FCC_norm = x;
x =  1 - abs(FD_FCC);   x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
FD_FCC_norm = x;
x =  ICCMed;    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
ICCMed_norm = x;
x =  ICCcontr;    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
ICCcontr_norm = x;

x = 1 - abs(xFDDVARS);    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
xFDDVARS_norm = x;
x =  1 - abs(xFD_FDDVARS);   x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
xFD_FDDVARS_norm = x;
x =  1 - abs(xFDFCmed);    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
xFDFCmed_norm = x;
x =  1 - abs(xFDFCM2);    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
xFDFCM2_norm = x;
x =  1 - abs(xFDFCdist);    x_raw = x(:,ind_raw);     x_raw_mean = mean(x_raw);    x_raw_std = std(x_raw);        x = x-x_raw_mean;    x = x./x_raw_std;
xFDFCdist_norm = x;

QC_signal = FCC_norm + ICCcontr_norm;  QC_signal = QC_signal./2;
QC_noise = FD_FCC_norm + xFDDVARS_norm +  xFD_FDDVARS_norm + xFDFCmed_norm + xFDFCM2_norm +  xFDFCdist_norm ;  QC_noise = QC_noise./6;
QC_all = QC_signal + QC_noise;  QC_all = QC_all./2;

N_PCAcomp_scale = [1:10,20:10:100,200:100:600]; N_PCAs = length(N_PCAcomp_scale);
color_k = [.1 .1 .1]; color_r = [.8 .0 .0];
color_b = [0.3   0.3    0.6]; color_g = [.0 0.9 0];

%% Plot scores for signal-related metrics


close all

if Section ==1
    figure('position', [  1282          83        1278        1273],  'Name','Scores for signal-related metrics','NumberTitle','off')
    for i = 1:4
        if i ==1
            x = FCC; label = 'FCC'; ylim_s =[46, 72];
        elseif i==2
            x = FD_FCC;  label = 'FD-FCC'; ylim_s = [-0.6  0.2];
        elseif i==3
            x = ICCMed;  label = 'MICC'; ylim_s =  [0.33, 0.85];
        elseif i == 4
            x = ICCcontr;  label =  'ICCC'; ylim_s =[18, 72];
        end
        
        x_ref = x(:,end-3:end);
        x(:, end-3:end) = [];
        
        % Without GSR   -----------------------------------
        subplot(5,2,(i-1)*2+1)
        for j = [1,3]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 1
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 3
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [1,2]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==1
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==2
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);         ylim(ylim_s)
        end
        if i==4, xlabel('# of PCA components'), end
        
        % With GSR   -----------------------------------
        subplot(5,2,(i-1)*2+2)
        for j = [2,4]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 2
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k);  hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 4
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r); hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [3,4]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==3
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==4
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);      ylim(ylim_s)
        end
        if i==4, xlabel('# of PCA components'), end
        
    end
    ax_c(1) = errorbar(-3,0,'LineWidth',2,'color',color_k); hold on
    ax_c(2) = errorbar(-3,0,'LineWidth',2,'color',color_r); hold on
    ax_c(3) = errorbar(-3,0,'LineWidth',2,'color',color_g); hold on
    ax_c(4) = errorbar(-3,0,'LineWidth',2,'color',color_b); hold on
    legend([ax_c],'Raw','FIX','WM','CSF')
end


%%  Plot scores for motion-related metrics


if Section == 1
    figure('position', [  1282          83        1278        1273],  'Name','Scores for motion-related metrics','NumberTitle','off')
    for i = 1:5
        if i ==1
            x = xFDDVARS; label = 'FDDVARS'; ylim_s =[-0.13, 0.42];
        elseif i==2
            x = xFD_FDDVARS;  label = 'FD-FDDVARS'; ylim_s =[-0.48 0.56];
        elseif i == 3
            x = xFDFCmed;  label = 'FDFC_{median}'; ylim_s =[0.11 0.27];
        elseif i == 4
            x = xFDFCM2;  label = 'FD-MFC'; ylim_s =[-0.17 0.6];
        else
            x = xFDFCdist;   label = 'FDFC_{dist}'; ylim_s =[-0.33 0.09];
        end
        
        
        x_ref = x(:,end-3:end);
        x(:, end-3:end) = [];
        
        % Without GSR   -----------------------------------
        subplot(5,2,(i-1)*2+1)
        for j = [1,3]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 1
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 3
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [1,2]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==1
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==2
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);          ylim(ylim_s)
        end
        if i==5, xlabel('# of PCA components'), end
        
        % With GSR   -----------------------------------
        subplot(5,2,(i-1)*2+2)
        for j = [2,4]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 2
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 4
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [3,4]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==3
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==4
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);       ylim(ylim_s)
        end
        if i==5, xlabel('# of PCA components'), end
        
    end
    ax_c(1) = errorbar(-3,0,'LineWidth',2,'color',color_k); hold on
    ax_c(2) = errorbar(-3,0,'LineWidth',2,'color',color_r); hold on
    ax_c(3) = errorbar(-3,0,'LineWidth',2,'color',color_g); hold on
    ax_c(4) = errorbar(-3,0,'LineWidth',2,'color',color_b); hold on
    legend([ax_c],'Raw','FIX','WM','CSF')
end


%%  Plot normalized scores for signal-related metrics

if Section ==1
    figure('position', [  1282          83        1278        1273],   'Name','Normalized scores for signal-related metrics','NumberTitle','off')
    for i = 1:4
        if i ==1
            x = FCC_norm; label = 'FCC'; ylim_s =[-1, 15];
        elseif i==2
            x = FD_FCC_norm;  label = 'FD-FCC'; ylim_s = [-1.35  3.3];
        elseif i==3
            x = ICCMed_norm;  label = 'MICC'; ylim_s =  [-12, 1.9];
        elseif i == 4
            x = ICCcontr_norm;  label =  'ICCC'; ylim_s =[-0.7, 5.2];
        end
        
        x_ref = x(:,end-3:end);
        x(:, end-3:end) = [];
        
        % Without GSR   -----------------------------------
        subplot(5,2,(i-1)*2+1)
        for j = [1,3]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 1
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 3
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [1,2]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==1
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==2
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);        ylim(ylim_s)
        end
        if i==4, xlabel('# of PCA components'), end
        
        % With GSR   -----------------------------------
        subplot(5,2,(i-1)*2+2)
        for j = [2,4]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 2
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 4
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [3,4]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==3
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==4
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);        ylim(ylim_s)
        end
        if i==4, xlabel('# of PCA components'), end
        
    end
    ax_c(1) = errorbar(-3,0,'LineWidth',2,'color',color_k); hold on
    ax_c(2) = errorbar(-3,0,'LineWidth',2,'color',color_r); hold on
    ax_c(3) = errorbar(-3,0,'LineWidth',2,'color',color_g); hold on
    ax_c(4) = errorbar(-3,0,'LineWidth',2,'color',color_b); hold on
    legend([ax_c],'Raw','FIX','WM','CSF')
end



%%  Plot normalized scores for motion-related metrics


if Section == 1
    figure('position', [  1282          83        1278        1273], 'Name', 'Normalized scores for motion-related metrics','NumberTitle','off')
    for i = 1:5
        if i ==1
            x = xFDDVARS_norm; label = 'FDDVARS'; ylim_s =[-1, 13];
        elseif i==2
            x = xFD_FDDVARS_norm;  label = 'FD-FDDVARS'; ylim_s =[-0.5 3.4];
        elseif i == 3
            x = xFDFCmed_norm;  label = 'FDFC_{median}'; ylim_s =[-2.7 1.3];
        elseif i == 4
            x = xFDFCM2_norm;  label = 'FD-MFC'; ylim_s =[-4.6 1.8];
        else
            x = xFDFCdist_norm;   label = 'FDFC_{dist}'; ylim_s =[-2.2 1.8];
        end
        
        
        x_ref = x(:,end-3:end);
        x(:, end-3:end) = [];
        
        % Without GSR   -----------------------------------
        subplot(5,2,(i-1)*2+1)
        for j = [1,3]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 1
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 3
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [1,2]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==1
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==2
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);          ylim(ylim_s)
        end
        if i==5, xlabel('# of PCA components'), end
        
        % With GSR   -----------------------------------
        subplot(5,2,(i-1)*2+2)
        for j = [2,4]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 2
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 4
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [3,4]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==3
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==4
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);      ylim(ylim_s)
        end
        if i==5, xlabel('# of PCA components'), end
    end
    ax_c(1) = errorbar(-3,0,'LineWidth',2,'color',color_k); hold on
    ax_c(2) = errorbar(-3,0,'LineWidth',2,'color',color_r); hold on
    ax_c(3) = errorbar(-3,0,'LineWidth',2,'color',color_g); hold on
    ax_c(4) = errorbar(-3,0,'LineWidth',2,'color',color_b); hold on
    legend([ax_c],'Raw','FIX','WM','CSF')
end



%%  Plot summarised scores


if Section == 1
    figure('position', [ 1385         368        1155         802 ], 'name', 'Summarized scores','NumberTitle','off')
    for i = 1:3
        if i ==1
            x = QC_signal; label = 'QC_{signal}'; ylim_s =[-1., 9.5];
        elseif i==2
            x = QC_noise;  label = 'QC_{motion}'; ylim_s =[-0.4 3.5];
        elseif i == 3
            x = QC_all;  label = 'CQC'; ylim_s =[-0.5 6.3];
        end
        
        x_ref = x(:,end-3:end);
        x(:, end-3:end) = [];
        
        % Without GSR   -----------------------------------
        subplot(3,2,(i-1)*2+1)
        for j = [1,3]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 1
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 3
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [1,2]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==1
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==2
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);           ylim(ylim_s)
        end
        if i==3, xlabel('# of PCA components'), end
        
        % With GSR   -----------------------------------
        subplot(3,2,(i-1)*2+2)
        for j = [2,4]
            ref_j = x_ref(:,j);        ref_j = repmat(ref_j,1,N_PCAs)    ;
            if j == 2
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_k), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_k); hold on
            elseif j == 4
                plot(N_PCAcomp_scale,mean(ref_j),'LineWidth',2, 'Color',color_r), hold on
                errorbar(500, mean(ref_j(:,1)),std(ref_j(:,1))/sqrt(10) ,'LineWidth',2,'color',color_r); hold on
            end
        end
        for b = [3,4]
            ind = (b-1)*N_PCAs + [1:N_PCAs];
            if b==3
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_g); hold on
            elseif b==4
                ax = errorbar(N_PCAcomp_scale, mean(x(:,ind)),std(x(:,ind))/sqrt(10) ,'LineWidth',2,'color',color_b); hold on
            end
            ax = gca;             set(ax,'Xscale','log');            ax.XGrid = 'on'; xlim([0 610]); ax.YGrid = 'on';
            ylabel(label);       ylim(ylim_s)
        end
        if i==3, xlabel('# of PCA components'), end
    end
    ax_c(1) = errorbar(-3,0,'LineWidth',2,'color',color_k); hold on
    ax_c(2) = errorbar(-3,0,'LineWidth',2,'color',color_r); hold on
    ax_c(3) = errorbar(-3,0,'LineWidth',2,'color',color_g); hold on
    ax_c(4) = errorbar(-3,0,'LineWidth',2,'color',color_b); hold on
    legend([ax_c],'Raw','FIX','WM','CSF')
end


%% S2:  Plot QC metrics for pipelines from literature

if Section ==2
    pos_a = [0.184589749577144,0.709264705882353,0.720410250422856,0.215735294117647];
    pos_b = [0.184589749577144,0.426942605937182,0.720410250422856,0.215735294117647];
    pos_c = [0.184589749577144,0.141957390146472,0.720410250422856,0.215735294117647];
    indPipel = [6:8,2,9:14,3,4,16:22,5];
    
    close all
    figure('position', [ 1505         499         343         751],  'Name','Evaluate pipelines used in the literature','NumberTitle','off')
    for i = 1:3
        if i ==1
            x = QC_signal; label = 'QC_{signal}'; ylim_s =[0, 10];
        elseif i==2
            x = QC_noise;  label = 'QC_{motion}'; ylim_s =[-0.5, 4];   % try with tranpose
        else
            x = QC_all;  label = 'CQC'; ylim_s =[0, 6.5];
        end
        
        x  = x(:,indPipel);
        
        ax(i) = subplot(3,1,i)  ;
        errorbar( mean(x),std(x)/sqrt(10) ,'o','LineWidth',1.5,'color','k','LineStyle','none','MarkerSize',2); hold on
        axg = gca;               axg.YGrid = 'on';
        ylabel(label);
        axg.XGrid = 'on';
        ylim(ylim_s)
        
    end
    linkaxes(ax,'x')
    xlim([0 23])
    xlabel('Pipeline')
    xlim([0 21])
    
    set(ax(1),'Position',pos_a)
    set(ax(2),'Position',pos_b)
    set(ax(3),'Position',pos_c)
    
end

%% S3:  Plot QC metrics for model-based NCTs

if Section == 3
    pos_1 = [0.130000000000000,0.819104097509062,0.775000000000000,0.157741935483871];
    pos_2 = [0.130000000000000,0.642437557543356,0.775000000000000,0.157741935483871];
    pos_3 = [0.130000000000000,0.460809402405487,0.775000000000000,0.157741935483871];
    pos_4 = [0.130000000000000,0.283435784851812,0.775000000000000,0.157741935483871];
    
    color_k = [.1 .1 .1];color_r = [.8 .0 .0];
    color_g = [0.1328    0.5430    0.1328];  color_b = [.0 .6 .9];
    
    close all
    load(['Struct_64.mat'])
    figure('position', [ 903   391   726   863],  'Name','Effect of model-based NCTs','NumberTitle','off')
    for i = 1:3
        if i ==1
            x = QC_signal; label = 'QC_{signal}'; ylim_s =[-1, 9.2];
        elseif i==2
            x = QC_noise;  label = 'QC_{motion}'; ylim_s =[-0.5, 3.6];   % try with tranpose
        elseif i == 3
            x = QC_all;  label = 'CQC'; ylim_s =[-0.5, 6.5];
        end
        
        ind =  ranges(1,1) : ranges(1,2);        x = x(:,ind);
        
        ax(i) = subplot(4,1,i)  ;
        errorbar( mean(x),std(x)/sqrt(10) ,'o','LineWidth',1.5,'color','k','LineStyle','none','MarkerSize',2); hold on
        axg = gca;               axg.YGrid = 'on';
        ylabel(label);     set(gca,'xtick',[])
        
        yl = [-10 20]; hold on;
        indModels = [8:8:64-8];
        nR = length(indModels);  cGrid = 0.5;
        for j = 1: nR
            xk = indModels(j);
            plot([xk+.5,xk+.5],[yl(1),yl(2)],'k-','LineWidth',0.01,'Color',[cGrid,cGrid,cGrid]);
        end
        ylim(ylim_s)
    end
    
    ax(4) = subplot(4,1,4);
    imagesc(struct_64')
    colormap_BW = [1 1 1 ; 0.5   0.5 0.5];         colormap(colormap_BW)
    axg =gca;
    
    set(gca,'xtick',[])
    yticks([1:6])
    yticklabels({'Cardiac ','Breathing','24 MPs','SLFOs','GS','WM^{200}'})
    set(gca,'TickLength',[0.0, 0.05])
    
    yl = [0.5 6.5]; hold on;
    indModels = [1:63];
    nR = length(indModels);  cGrid = 0.7;
    for j = 1: nR
        xk = indModels(j);
        plot([xk+.5,xk+.5],[yl(1),yl(2)],'k-','LineWidth',0.01,'Color',[cGrid,cGrid,cGrid]);
    end
    
    indModels = [8:8:64-8];
    nR = length(indModels);  cGrid = 0.3;
    for j = 1: nR
        xk = indModels(j);
        plot([xk+.5,xk+.5],[yl(1),yl(2)],'k-','LineWidth',0.01,'Color',[cGrid,cGrid,cGrid]);
    end
    ylim(yl)
    
    indModels = [1:5];
    nR = length(indModels);  cGrid = 0.3;
    for j = 1: nR
        yk = indModels(j);
        plot( [0,65], [yk+.5,yk+.5],'k-','LineWidth',0.01,'Color',[cGrid,cGrid,cGrid]);
    end
    ylim(yl)
    
    xlabel('Model')
    linkaxes(ax,'x')
    xlim([0.5 64.5])
    
    set(ax(1),'Position',pos_1)
    set(ax(2),'Position',pos_2)
    set(ax(3),'Position',pos_3)
    set(ax(4),'Position',pos_4)    
    
end

%% S4: Plots for censoring

if Section == 4
    FD_thresh_list = [0.1:0.05:1.0,inf];  % FD
    DVARS_thresh_list = [0.5:0.5:20,inf];        % DVARS
    indFD = [20, 19, 15, 9, 5,4,3,2];  FD_thresh_list(indFD);
    indDVARS = [41, 40, 20, 10,  4,3,2,1];  DVARS_thresh_list(indDVARS);
    
    
    close all
    figure('position', [   1436         678        1017         562 ], 'Name','Effect of scrubbing','NumberTitle','off')
    for i = 1:3
        if i ==1
            x = QC_signal; label = 'QC_{signal}'; ylim_s =[5.5, 9.5];
        elseif i==2
            x = QC_noise;  label = 'QC_{motion}'; ylim_s =[0.5, 3.5];
        elseif i == 3
            x = QC_all;  label = 'CQC'; ylim_s =[3.2, 6.2];
        end
        
        ind =  ranges(1,1) : ranges(1,2);        xFD = x(:,ind(indFD));       % FD
        ind =  ranges(2,1) : ranges(2,2);        xDVARS = x(:,ind(indDVARS));     % DVARS
        
        ax1(i) = subplot(3,2,(i-1)*2+1)  ;
        errorbar( mean(xFD),std(xFD)/sqrt(10) ,'o','LineWidth',1.5,'color','k','LineStyle','none','MarkerSize',2); hold on
        ax1(i).YGrid = 'on';
        
        ylabel(label);       ylim(ylim_s)
        xticks([1:length(indFD)]) ,       xticklabels({'+\infty','1.00','0.80','0.50','0.30','0.25','0.20','0.15'})
        if i == 3,     xlabel('FD_{thr} (mm)'),    end
        
        ax2(i) = subplot(3,2,(i-1)*2+2);
        errorbar( mean(xDVARS),std(xDVARS)/sqrt(10) ,'o','LineWidth',1.5,'color','k','LineStyle','none','MarkerSize',2); hold on
        ax2(i).YGrid = 'on';
        ylabel(label);       ylim(ylim_s)
        xticks([1:length(indDVARS)]) ,       xticklabels({'+\infty','20','10','5','2','1.5','1','0.5'})
        
        if i == 3,         xlabel('DVARS_{thr} (MAD)'),      end
    end
    linkaxes(ax1,'x'),         subplot(ax1(1)),     xlim([0.5 8.5])
    linkaxes(ax2,'x'),   subplot(ax2(1)),  xlim([0.5 8.5])
end

%%  S5: Plot QC metrics for the effect of low-pass filtering

if Section ==5
    close all
    figure('position', [ 1468         634         508         557  ], 'Name','Effect of low-pass filtering','NumberTitle','off' )
    for i = 1:3
        if i == 1
            x = QC_signal; label = 'QC_{signal}'; ylim_s =[6.7, 9.5];
        elseif i==2
            x = QC_noise;  label = 'QC_{motion}'; ylim_s =[2.3, 3.5];   % try with tranpose
        elseif i == 3
            x = QC_all;  label = 'CQC'; ylim_s =[4.8, 6.3];
        end
        
        ind =  ranges(1,1) : ranges(2,2);        x = x(:,ind);   x1 = x(:,1:10); x2 = x(:,11:20);
        ax(i) = subplot(3,1,i)  ;
        errorbar( mean(x2),std(x2)/sqrt(10) ,'o','LineWidth',1.5,'color','k','LineStyle','none','MarkerSize',2); hold on
        axg = gca;               axg.YGrid = 'on';
        ylabel(label);     set(gca,'xtick',[]) ,           ylim(ylim_s)
        
        set(gca,'xtick',[]),         xticks([1:10])
        xticklabels({'+\infty','0.6','0.5','0.4','0.3','0.2','0.1','0.08','0.05','0.01'})
        if i ==3,     xlabel('Cut-off frequency (Hz)')   , end
        
    end
    
    linkaxes(ax,'x')
    xlim([0.5 9.5])
end


%%  ==============================












