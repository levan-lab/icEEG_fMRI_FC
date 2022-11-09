% computes dynamic fc for fMRI and icEEG and stores the values in a variable for R
clear all
sub = 'sub-049';
run = '02';

projdir = '/Users/kristinasabaroedin/Documents/icEEG_fMRI/';

% we want to lag BOLD fc by 5 seconds
% because TR = 1.5, 5 seconds do not align to a start of a volume
% so I shift BOLD by six seconds, and shift EEG by 1 sec

% index for BOLD window
bold_dyn_secs = linspace(60,1200,20);
bold_dyn_secs = bold_dyn_secs + 6;
bold_dyn_secs = [6, bold_dyn_secs];
bold_dyn_secs = bold_dyn_secs(1,1:20);
bold_dyn_vols = bold_dyn_secs./1.5;


% index for EEG window
a = linspace(61,1201,20);
time_intervals = [1, a];
time_intervals = time_intervals(1,1:20);
clear a

% initiate structure for storing conn metrics
conn_metrics_120s = {};


% compute dynamic fc
cd([projdir,sub,'/conn_analysis/'])

mkdir('sliding_win')
mkdir('sliding_win/120s')

load(['run',run,'_24hmp2p_butter_rois_fc_select.mat'])
numRois = length(roi_pairs);

tscorr = {};


for i = 1:18
    j = i+2;
    tscorr_fcmatrix{i} = corr(ts(bold_dyn_vols(i):bold_dyn_vols(j),:));
end


bold_fc = [];
    
for j = 1:18
    corrt = tscorr_fcmatrix{j}; 
    B = ones(numRois,numRois);
    tri = triu(corrt);
    bold_fc(:,j) = tri(logical(triu(B,1)));
    
end

conn_metrics_120s.boldfcmatrix = tscorr_fcmatrix;


conn_metrics_120s.boldfc = bold_fc;
conn_metrics_120s.std_boldfc = std(conn_metrics_120s.boldfc')';
conn_metrics_120s.bold_dyn_secs = bold_dyn_secs;
conn_metrics_120s.bold_dyn_vols = bold_dyn_vols;

clear tscorr ts l j i a bold_fc contacts_corr corrt


%%%%%%%%%%%

cd([projdir,sub,'/icEEG/'])
load([sub,'_run',run,'_icEEG_filt.mat'])

time = {};
time{1,1} = linspace(1,1,2400000);

data = {};
data.time = time;
data.fsample = fs;
data.label = label;
data.trial{1,1} = eeg_filt;



%% 
cfg = [];
cfg.length = 1;
cfg.overlap = 0;
data_segmented = ft_redefinetrial(cfg, data)
%% 

freq = {};
fd = {};

for i = 1:18
    j = i+2;
    cfg = [];
    cfg.output = 'powandcsd'
    cfg.method     = 'mtmfft';
    cfg.foilim     = [1 100];
    cfg.taper       = 'dpss'; 
    cfg.tapsmofrq = 5;
    cfg.trials = time_intervals(i):time_intervals(j)
    cfg.keeptrials = 'yes';
    freq{i}  = ft_freqanalysis(cfg, data_segmented)
    %%%
    cfg            = [];
    cfg.method     = 'coh';
    fd{i} = ft_connectivityanalysis(cfg, freq{i});
    %%%
    freq{i}.powspctrm = permute(freq{i}.powspctrm, [2,3,1]);
    freq{i}.dimord = 'chan_freq_time';
    
end

for j = 1:18
    for i = 1:length(fd{1}.labelcmb)

        delta(i,:) = fd{j}.cohspctrm(i,1:3);
        delta_mean(i,j) = mean(delta(i,:));

        theta(i,:) = fd{j}.cohspctrm(i,4:7);
        theta_mean(i,j) = mean(theta(i,:));

        alpha(i,:) = fd{j}.cohspctrm(i,8:12);
        alpha_mean(i,j) = mean(alpha(i,:));

        beta(i,:) = fd{j}.cohspctrm(i,13:30);
        beta_mean(i,j) = mean(beta(i,:));

        gamma_low(i,:) = fd{j}.cohspctrm(i,31:50);
        gamma_low_mean(i,j) = mean(gamma_low(i,:));

        gamma_high(i,:) = fd{j}.cohspctrm(i,51:100);
        gamma_high_mean(i,j) = mean(gamma_high(i,:));

        broadband(i,:) = fd{j}.cohspctrm(i,1:100);
        broadband_mean(i,j) = mean(broadband(i,:));

    end
end

load([projdir,sub,'/conn_analysis/',sub,'_run',run,'_static_conn.mat'])


conn_metrics_120s.freq = freq;
conn_metrics_120s.fd = fd;
conn_metrics_120s.delta_coh = delta_mean;
conn_metrics_120s.theta_coh = theta_mean;
conn_metrics_120s.alpha_coh = alpha_mean;
conn_metrics_120s.beta_coh = beta_mean;
conn_metrics_120s.gammalow_coh = gamma_low_mean;
conn_metrics_120s.gammahigh_coh = gamma_high_mean;
conn_metrics_120s.broadband_coh = broadband_mean;
conn_metrics_120s.labelcmb = fd{1}.labelcmb;
conn_metrics_120s.seconds = time_intervals(1,3:end);


%% get conn n x n matrix for coherence double

a(1,1) = 1;
b(1,1) = numRois-1;

for i = 2:numRois
    
   a(i,1) = (a(i-1,1)+1) + (numRois-i);
   b(i,1) = a(i,1) + (numRois-(i+1));
  
end

coh_matrix_delta = eye(numRois);
coh_matrix_theta = eye(numRois);
coh_matrix_alpha = eye(numRois);
coh_matrix_beta = eye(numRois);
coh_matrix_gammalow = eye(numRois);
coh_matrix_gammahigh = eye(numRois);
coh_matrix_broadband = eye(numRois);

for j = 1:18
    coh =  conn_metrics_120s.delta_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_delta(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_delta(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    delta_coh_matrix{j} = coh_matrix_delta;
end


for j = 1:18
    coh =  conn_metrics_120s.theta_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_theta(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_theta(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    theta_coh_matrix{j} = coh_matrix_theta;
end

for j = 1:18
    coh =  conn_metrics_120s.alpha_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_alpha(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_alpha(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    alpha_coh_matrix{j} = coh_matrix_alpha;
end

for j = 1:18
    coh =  conn_metrics_120s.beta_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_beta(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_beta(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    beta_coh_matrix{j} = coh_matrix_beta;
end


for j = 1:18
    coh =  conn_metrics_120s.gammalow_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_gammalow(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_gammalow(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    gammalow_coh_matrix{j} = coh_matrix_gammalow;
end

for j = 1:18
    coh =  conn_metrics_120s.gammahigh_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_gammahigh(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_gammahigh(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    gammahigh_coh_matrix{j} = coh_matrix_gammahigh;
end

for j = 1:18
    coh =  conn_metrics_120s.broadband_coh(:,j);
    for g = 1:numRois-1
        coh_matrix_broadband(g,g+1:end) = coh(a(g,1):b(g,1));
        coh_matrix_broadband(g+1:end,g) = coh(a(g,1):b(g,1));
    end
    
    broadband_coh_matrix{j} = coh_matrix_broadband;
end


conn_metrics_120s.coh_matrix_delta = delta_coh_matrix;
conn_metrics_120s.coh_matrix_theta = theta_coh_matrix;
conn_metrics_120s.coh_matrix_alpha = alpha_coh_matrix;
conn_metrics_120s.coh_matrix_beta = beta_coh_matrix;
conn_metrics_120s.coh_matrix_gammalow = gammalow_coh_matrix;
conn_metrics_120s.coh_matrix_gammahigh = gammahigh_coh_matrix;
conn_metrics_120s.coh_matrix_broadband = broadband_coh_matrix;

clear coh_matrix_delta coh_matrix_theta coh_matrix_alpha 
clear coh_matrix_beta coh_matrix_gammalow coh_matrix_gammahigh coh_matrix_broadband 
clear delta_coh_matrix theta_coh_matrix alpha_coh_matrix beta_coh_matrix gammalow_coh_matrix
clear gammahigh_coh_matrix broadband_coh_matrix

%%%%%

std_delta_coh = std(conn_metrics_120s.delta_coh')';
std_theta_coh = std(conn_metrics_120s.theta_coh')';
std_alpha_coh = std(conn_metrics_120s.alpha_coh')';
std_beta_coh = std(conn_metrics_120s.beta_coh')';
std_gammalow_coh = std(conn_metrics_120s.gammalow_coh')';
std_gammahigh_coh = std(conn_metrics_120s.gammahigh_coh')';
std_broadband_coh = std(conn_metrics_120s.broadband_coh')';

conn_metrics_120s.std_coh = table(std_delta_coh, std_theta_coh, std_alpha_coh, std_beta_coh,...
    std_gammalow_coh, std_gammahigh_coh, std_broadband_coh);



clear std_delta_coh std_theta_coh std_alpha_coh std_beta_coh
clear std_gammalow_coh std_gammahigh_coh std_broadband_coh

% get correlation between bold fc and coherence
for i = 1:length(conn_metrics_120s.labelcmb)
    [r_bold_delta(i,:), p_bold_delta(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.delta_coh(i,:))');
    [r_bold_theta(i,:), p_bold_theta(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.theta_coh(i,:))');
    [r_bold_alpha(i,:), p_bold_alpha(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.alpha_coh(i,:))');
    [r_bold_beta(i,:), p_bold_beta(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.beta_coh(i,:))');
    [r_bold_gammalow(i,:), p_bold_gammalow(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.gammalow_coh(i,:))');
    [r_bold_gammahigh(i,:), p_bold_gammahigh(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.gammahigh_coh(i,:))');
    [r_bold_broadband(i,:), p_bold_broadband(i,:)] = corr(atanh(conn_metrics_120s.boldfc(i,:))',atanh(conn_metrics_120s.broadband_coh(i,:))');
end



 conn_metrics_120s.bold_coh_corr = table(r_bold_delta, p_bold_delta, r_bold_theta, p_bold_theta, r_bold_alpha,  p_bold_alpha,...
    r_bold_beta, p_bold_beta, r_bold_gammalow, p_bold_gammalow, r_bold_gammahigh, p_bold_gammahigh,...
    r_bold_broadband, p_bold_broadband);

clear r_bold_delta p_bold_delta r_bold_theta  p_bold_theta r_bold_alpha p_bold_alpha
clear r_bold_beta p_bold_beta r_bold_gammalow p_bold_gammalow r_bold_gammahigh  
clear p_bold_gammahigh r_bold_broadband p_bold_broadband


%% 
%%% compute power spectrum and calculate correlation of power between
%%% contacts

delta = {};
theta = {};
alpha = {};
beta = {};
gammalow = {};
gammahigh = {};
broadband = {};

for i = 1:18
    for k = 1:120 % seconds
        delta{i}(:,k) = mean(freq{1,i}.powspctrm(:,1:3,k),2);
        theta{i}(:,k) = mean(freq{1,i}.powspctrm(:,4:7,k),2);
        alpha{i}(:,k) = mean(freq{1,i}.powspctrm(:,8:12,k),2);
        beta{i}(:,k) = mean(freq{1,i}.powspctrm(:,13:30,k),2);
        gammalow{i}(:,k) = mean(freq{1,i}.powspctrm(:,31:50,k),2);
        gammahigh{i}(:,k) = mean(freq{1,i}.powspctrm(:,51:100,k),2);
        broadband{i}(:,k) = mean(freq{1,i}.powspctrm(:,1:100,k),2);
    end
end



delta_pow_corr_matrix = {};
theta_pow_corr_matrix = {};
alpha_pow_corr_matrix = {};
beta_pow_corr_matrix = {};
gammalow_pow_corr_matrix = {};
gammahigh_pow_corr_matrix = {};
broadband_pow_corr_matrix = {};

% compute correlations for power between all contacts
for i = 1:18
    delta_pow_corr_matrix{i} = corr(delta{1,i}');
    theta_pow_corr_matrix{i} = corr(theta{1,i}');
    alpha_pow_corr_matrix{i} = corr(alpha{1,i}');
    beta_pow_corr_matrix{i} = corr(beta{1,i}');
    gammalow_pow_corr_matrix{i} = corr(gammalow{1,i}');
    gammahigh_pow_corr_matrix{i} = corr(gammahigh{1,i}');
    broadband_pow_corr_matrix{i} = corr(broadband{1,i}');
end


% now we want to change the power conn matrix to a double

delta_pow_corr = [];
theta_pow_corr = [];
alpha_pow_corr = [];
beta_pow_corr = [];
gammalow_pow_corr = [];
gammahigh_pow_corr = [];
broadband_pow_corr = [];


for j = 1:18
    corrpow =  delta_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    delta_pow_corr(:,j) = tri(logical(triu(B,1)));
end


for j = 1:18
    corrpow =  theta_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    theta_pow_corr(:,j) = tri(logical(triu(B,1)));  
end

for j = 1:18
    corrpow =  alpha_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    alpha_pow_corr(:,j) = tri(logical(triu(B,1)));
end


for j = 1:18
    corrpow =  beta_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    beta_pow_corr(:,j) = tri(logical(triu(B,1)));
end

for j = 1:18
    corrpow =  gammalow_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    gammalow_pow_corr(:,j) = tri(logical(triu(B,1)));
    
end


for j = 1:18
    corrpow =  gammahigh_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    gammahigh_pow_corr(:,j) = tri(logical(triu(B,1)));
  
end


for j = 1:18
    corrpow =  broadband_pow_corr_matrix{j};
    B = ones(numRois,numRois);
    tri = triu(corrpow);
    broadband_pow_corr(:,j) = tri(logical(triu(B,1)));

end

clear d corrpow


conn_metrics_120s.powspec_delta = delta;
conn_metrics_120s.powspec_theta = theta;
conn_metrics_120s.powspec_alpha = alpha;
conn_metrics_120s.powspec_beta = beta;
conn_metrics_120s.powspec_gammalow = gammalow;
conn_metrics_120s.powspec_gammahigh = gammahigh;
conn_metrics_120s.powspec_broadband = broadband;

conn_metrics_120s.powspec_corr_delta = delta_pow_corr;
conn_metrics_120s.powspec_corr_theta = theta_pow_corr;
conn_metrics_120s.powspec_corr_alpha = alpha_pow_corr;
conn_metrics_120s.powspec_corr_beta = beta_pow_corr;
conn_metrics_120s.powspec_corr_gammalow = gammalow_pow_corr;
conn_metrics_120s.powspec_corr_gammahigh = gammahigh_pow_corr;
conn_metrics_120s.powspec_corr_broadband = broadband_pow_corr;

conn_metrics_120s.powspec_corr_matrix_delta = delta_pow_corr_matrix;
conn_metrics_120s.powspec_corr_matrix_theta = theta_pow_corr_matrix;
conn_metrics_120s.powspec_corr_matrix_alpha = alpha_pow_corr_matrix;
conn_metrics_120s.powspec_corr_matrix_beta = beta_pow_corr_matrix;
conn_metrics_120s.powspec_corr_matrix_gammalow = gammalow_pow_corr_matrix;
conn_metrics_120s.powspec_corr_matrix_gammahigh = gammahigh_pow_corr_matrix;
conn_metrics_120s.powspec_corr_matrix_broadband = broadband_pow_corr_matrix;

clear delta theta alpha beta gammalow gammahigh broadband
clear delta_corr theta_corr alpha_corr beta_corr gammalow_corr gammahigh_corr broadband_corr
clear delta_pow_corr theta_pow_corr alpha_pow_corr beta_pow_corr gammalow_pow_corr gammahigh_pow_corr broadband_pow_corr
clear delta_pow_corr_matrix theta_pow_corr_matrix alpha_pow_corr_matrix beta_pow_corr_matrix gammalow_pow_corr_matrix
clear gammahigh_pow_corr_matrix broadband_pow_corr_matrix



save(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',sub,'/conn_analysis/sliding_win/120s/',sub,'_run',run,'_conn_metrics_120s.mat'],'-v7.3')

% ft_singleplotER(conn_metrics_120s.freq{1,1}.cfg, conn_metrics_120s.freq{1,1})
cd(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',sub,'/conn_analysis/sliding_win/120s/'])
for i = 1:18
    ft_singleplotER(conn_metrics_120s.freq{1,1}.cfg, conn_metrics_120s.freq{1,i})
    print -dpng
end