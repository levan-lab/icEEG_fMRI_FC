% script for static spectral analysis

clear all
addpath('/Users/kristinasabaroedin/Applications/fieldtrip-20210713')
ft_defaults

subject = 'sub-049';
run = '02';

% first, segment continuous data so fieltrip understands input
cd(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject,'/icEEG/'])
load(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject,'/icEEG/',subject,'_run',run,'_icEEG_filt.mat'])

time = {};
time{1,1} = linspace(1,1,2400000);

data = {};
data.time = time;
data.fsample = fs;
data.label = label;
data.trial{1,1} = eeg_filt;


%% 
% cfg = [];
% cfg.viewmode = 'vertical'
% ft_databrowser(cfg, data)
%% 

%% 
cfg = [];
cfg.length = 1.5;
cfg.overlap = 0;
data_segmented = ft_redefinetrial(cfg, data)
%% 


%%%%%%

cfg = [];
cfg.output = 'powandcsd'
cfg.method     = 'mtmfft';
cfg.foilim     = [1 100];
cfg.taper       = 'dpss'; 
cfg.tapsmofrq = 5;
cfg.keeptrials = 'yes';
freq_segmented = ft_freqanalysis(cfg, data_segmented)

begsample = data_segmented.sampleinfo(:,1);
endsample = data_segmented.sampleinfo(:,2);
time = ((begsample+endsample)/2) / data_segmented.fsample;

freq_continuous           = freq_segmented;
freq_continuous.powspctrm = permute(freq_segmented.powspctrm, [2, 3, 1]);
freq_continuous.dimord    = 'chan_freq_time'; % it used to be 'rpt_chan_freq'
freq_continuous.time      = time;             % add the description of the time dimension



%% 
cfg            = [];
cfg.method     = 'coh';
fd             = ft_connectivityanalysis(cfg, freq_segmented);
%% 



numRois = length(label); %number of channels
comb = numRois*(numRois-1)/2; % all combination of pairs



for i = 1:length(fd.labelcmb)
    delta(i,:) = fd.cohspctrm(i,1:3);
    delta_mean(i,:) = mean(delta(i,:));
    
    theta(i,:) = fd.cohspctrm(i,4:7);
    theta_mean(i,:) = mean(theta(i,:));
    
    alpha(i,:) = fd.cohspctrm(i,8:12);
    alpha_mean(i,:) = mean(alpha(i,:));
    
    beta(i,:) = fd.cohspctrm(i,13:30);
    beta_mean(i,:) = mean(beta(i,:));
    
    gamma_low(i,:) = fd.cohspctrm(i,31:50);
    gamma_low_mean(i,:) = mean(gamma_low(i,:));
    
    gamma_high(i,:) = fd.cohspctrm(i,51:100);
    gamma_high_mean(i,:) = mean(gamma_high(i,:));
    
    broadband(i,:) = fd.cohspctrm(i,1:100);
    broadband_mean(i,:) = mean(broadband(i,:));
 
end

coh_pow_freq = {};
coh_pow_freq.delta_mean = delta_mean;
coh_pow_freq.theta_mean = theta_mean;
coh_pow_freq.alpha_mean = alpha_mean;
coh_pow_freq.beta_mean = beta_mean;
coh_pow_freq.gammalow_mean = gamma_low_mean;
coh_pow_freq.gammahigh_mean = gamma_high_mean;
coh_pow_freq.broadband_mean = broadband_mean;

clear delta delta_mean theta theta_mean alpha alpha_mean beta beta_mean gamma_low gamma_low_mean gamma_high gamma_high_mean broadband broadband_mean



a(1,1) = 1;
b(1,1) = numRois-1;

for i = 2:numRois
    
   a(i,1) = (a(i-1,1)+1) + (numRois-i);
   b(i,1) = a(i,1) + (numRois-(i+1));
  
end

coh_conn_matrix = {};
coh_conn_matrix.delta = eye(numRois);
coh_conn_matrix.theta = eye(numRois);
coh_conn_matrix.alpha = eye(numRois);
coh_conn_matrix.beta = eye(numRois);
coh_conn_matrix.gammalow = eye(numRois);
coh_conn_matrix.gammahigh = eye(numRois);
coh_conn_matrix.broadband = eye(numRois);



for g = 1:numRois-1 % number of electrodes - 1
    coh_conn_matrix.delta(g,g+1:end) = coh_pow_freq.delta_mean(a(g,1):b(g,1));
    coh_conn_matrix.delta(g+1:end,g) = coh_pow_freq.delta_mean(a(g,1):b(g,1));

    coh_conn_matrix.theta(g,g+1:end) = coh_pow_freq.theta_mean(a(g,1):b(g,1));
    coh_conn_matrix.theta(g+1:end,g) = coh_pow_freq.theta_mean(a(g,1):b(g,1));

    coh_conn_matrix.alpha(g,g+1:end) = coh_pow_freq.alpha_mean(a(g,1):b(g,1));
    coh_conn_matrix.alpha(g+1:end,g) = coh_pow_freq.alpha_mean(a(g,1):b(g,1));

    coh_conn_matrix.beta(g,g+1:end) = coh_pow_freq.beta_mean(a(g,1):b(g,1));
    coh_conn_matrix.beta(g+1:end,g) = coh_pow_freq.beta_mean(a(g,1):b(g,1));

    coh_conn_matrix.gammalow(g,g+1:end) = coh_pow_freq.gammalow_mean(a(g,1):b(g,1));
    coh_conn_matrix.gammalow(g+1:end,g) = coh_pow_freq.gammalow_mean(a(g,1):b(g,1));

    coh_conn_matrix.gammahigh(g,g+1:end) = coh_pow_freq.gammahigh_mean(a(g,1):b(g,1));
    coh_conn_matrix.gammahigh(g+1:end,g) = coh_pow_freq.gammahigh_mean(a(g,1):b(g,1));

    coh_conn_matrix.broadband(g,g+1:end) = coh_pow_freq.broadband_mean(a(g,1):b(g,1));
    coh_conn_matrix.broadband(g+1:end,g) = coh_pow_freq.broadband_mean(a(g,1):b(g,1));

end

save(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject,'/conn_analysis/',subject,'_run',run,'_filteeg_powandcsd_coh_multitaper5Hz.mat'],'fd','freq_segmented','freq_continuous','coh_pow_freq','coh_conn_matrix','label')