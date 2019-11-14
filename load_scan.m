
function [ROI_data, GMm, WMm, WMpca, CSFm, CSFpca, FD,movRegr, RETR_RespRegr, RETR_CardRegr,GBF,DVARS,PCAexpl, GSstats] = load_scan(subject,task, flag_physio)


baseDir='E:\CloudStation\HCP\RawData/';
filepath_movRegr=[baseDir,'/Physio/',subject,'_',task,'/Movement_Regressors_dt.txt'];
filepath_MRacq=[baseDir,'/Physio/',subject,'_',task,'/phys.mat'];  load(filepath_MRacq,'Fs','trig','TR'); Ts = 1/Fs;
load([baseDir,'/Physio/',subject,'_',task,'/Phys_sum.mat']);

filepath_input=[baseDir,'Atlas/',subject,'_',task,'/'];
load([filepath_input,'TissueBasedRegressors_1199.mat'])

load( 'Atlases/Gordon/Gordon333.mat')
load([filepath_input,'ROI_data_Gordon_333_surf.mat']),  ROI_data = ROI_data(:, parcels_rem+1:end); NC = size(ROI_data,2);

%%  Create timeline   ---------------

volDel=40;
ind_BOLD=find(trig==1);     trig(ind_BOLD(1:volDel))=0;    ind_BOLD=find(trig==1);
time = 0:Ts:(length(trig)-1)*Ts;
timeMR=time(find(trig==1));

ind_BOLD_10=zeros(length(timeMR),1);
for i=1:length(timeMR)
    [~,ind]=min(abs(time_10-timeMR(i)));   %% consider number of volumes before task
    ind_BOLD_10(i)=ind;
end

ROI_data(1:volDel,:) = [];  [NV,NC]=size(ROI_data);
DVARS = DVARS(1+volDel:end);

if ~(NV==sum(trig)), disp('Check NV!!!'), end

%   Load motion parameters (calculate FD)   ---------------------

movRegr=load(filepath_movRegr); movRegr=(movRegr(1+volDel:end,:));
movRegr=[movRegr, movRegr.^2]; NR_motion=size(movRegr,2);
Mov = abs(movRegr(:,7:12));    Mov1 = sum(Mov(:,1:3),2) ;    Mov2 =  50 *  sum(Mov(:,4:6),2) * (2*pi()/360) ;    FD = Mov1 + Mov2;  movRegr = zscore(movRegr);

% fprintf('Mean FD: %3.2f \n', mean(FD))

%% Preprocess parcel timeseries and nuisance regressors

HPF_f = 0.008;
[filt_b,filt_a] = butter(2,HPF_f*2*TR,'high');

ROI_data = filtfilt(filt_b,filt_a,ROI_data);
movRegr = filtfilt(filt_b,filt_a,movRegr ) ;


PCAind = 1:size(WB.PCAcomp,2);
TBR_list = {'GM','WM','CSF'};
ind = 1+volDel:NV+volDel;
GSstats.mean = mean(GM.MA(ind));
GSstats.std = std(GM.MA(ind));

for i=1:length(TBR_list)
    eval(sprintf('x = %s ;', TBR_list{i}));
    MA = x.MA(ind);
    x.MA =  filtfilt(filt_b,filt_a,MA );
    
    x.PCAcomp=x.PCAcomp(ind,PCAind);
    x.PCAcomp = filtfilt(filt_b,filt_a,x.PCAcomp ) ;
    eval(sprintf('%s = x ;', TBR_list{i}));
end

%%  Canonical correlation bw regressors   ----------------------------------

movRegr = detrend(movRegr,'linear');
GMm = detrend(GM.MA,'linear');
WMm = detrend(WM.MA,'linear');
WMpca = detrend(WM.PCAcomp,'linear');
CSFm = detrend(CSF.MA,'linear');
CSFpca = detrend(CSF.PCAcomp,'linear');

PCAexpl = [WM.explained, CSF.explained];

%%  -------------------------

RETR_RespRegr = [];  RETR_CardRegr = []; GBF= [];

if flag_physio == 1
    %     RETR_RespRegr=RETR_Resp_regressors(resp,2,Fs); RETR_RespRegr = RETR_RespRegr(ind_BOLD,:);
    %     RETR_CardRegr=RETR_Card_regressors_v2(time,PPGlocs,3);   RETR_CardRegr = RETR_CardRegr(ind_BOLD,:);
    
    load([baseDir,'/Physio/',subject,'_',task,'/RETROICOR.mat']);
    RETR_RespRegr = RETR_RespRegr(volDel+1:end,1:6);
    RETR_CardRegr = RETR_CardRegr(volDel+1:end,:);
    
    load([baseDir,'/Physio/',subject,'_',task,'/Global_blood_flow.mat']);
    
    RETR_RespRegr = filtfilt(filt_b,filt_a,RETR_RespRegr); RETR_RespRegr = detrend(RETR_RespRegr,'linear');
    RETR_CardRegr = filtfilt(filt_b,filt_a,RETR_CardRegr); RETR_CardRegr = detrend(RETR_CardRegr,'linear');
    GBF = [GBF_HR,GBF_RF];     GBF = filtfilt(filt_b,filt_a,GBF);         GBF = detrend(GBF,'linear');
    
end

end

