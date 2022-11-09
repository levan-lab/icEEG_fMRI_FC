%%% compute power spectrum and calculate correlation of power between contacts

clear
subject = 'sub-049';
run = '02';

cd(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject,'/conn_analysis/']);
load([subject,'_run',run,'_filtEEG_powandcsd_coh_multitaper5Hz.mat'],'freq_continuous','label')
numRois = length(label)
%% 


delta = [];
theta = [];
alpha = [];
beta = [];
gammalow = [];
gammahigh = [];
broadband = [];

for i = 1:length(freq_continuous.label)
    a = mean(freq_continuous.powspctrm(i,1:3,:));
    delta(:,i) = reshape(a,1200,1);
    clear a
    a = mean(freq_continuous.powspctrm(i,4:7,:));
    theta(:,i) = reshape(a,1200,1);
    clear a
    a = mean(freq_continuous.powspctrm(i,8:12,:));
    alpha(:,i) = reshape(a,1200,1);
    clear a
    a = mean(freq_continuous.powspctrm(i,13:30,:));
    beta(:,i) = reshape(a,1200,1);
    clear a
    a = mean(freq_continuous.powspctrm(i,31:50,:));
    gammalow(:,i) = reshape(a,1200,1);
    clear a
    a = mean(freq_continuous.powspctrm(i,51:100,:));
    gammahigh(:,i) = reshape(a,1200,1);
    clear a
    a = mean(freq_continuous.powspctrm(i,1:100,:));
    broadband(:,i) = reshape(a,1200,1);
end

clear a

delta_pow_corr_matrix = corr(delta);
theta_pow_corr_matrix = corr(theta);
alpha_pow_corr_matrix = corr(alpha);
beta_pow_corr_matrix = corr(beta);
gammalow_pow_corr_matrix = corr(gammalow);
gammahigh_pow_corr_matrix = corr(gammahigh);
broadband_pow_corr_matrix = corr(broadband);


% now we want to change the power conn matrix to a double

delta_pow_corr = [];
theta_pow_corr = [];
alpha_pow_corr = [];
beta_pow_corr = [];
gammalow_pow_corr = [];
gammahigh_pow_corr = [];
broadband_pow_corr = [];



B = ones(numRois,numRois);

tri = triu(delta_pow_corr_matrix);
delta_pow_corr = tri(logical(triu(B,1)));  
clear tri

tri = triu(theta_pow_corr_matrix);
theta_pow_corr = tri(logical(triu(B,1)));  
clear tri

tri = triu(alpha_pow_corr_matrix);
alpha_pow_corr = tri(logical(triu(B,1)));  
clear tri

tri = triu(beta_pow_corr_matrix);
beta_pow_corr = tri(logical(triu(B,1)));  
clear tri

tri = triu(gammalow_pow_corr_matrix);
gammalow_pow_corr = tri(logical(triu(B,1)));  
clear tri

tri = triu(gammahigh_pow_corr_matrix);
gammahigh_pow_corr = tri(logical(triu(B,1)));  
clear tri


tri = triu(broadband_pow_corr_matrix);
broadband_pow_corr = tri(logical(triu(B,1)));  
clear tri




pow_spec_mean = {};
pow_spec_mean.delta = delta;
pow_spec_mean.theta = theta;
pow_spec_mean.alpha = alpha;
pow_spec_mean.beta = beta;
pow_spec_mean.gammalow = gammalow;
pow_spec_mean.gammahigh = gammahigh;
pow_spec_mean.broadband = broadband;

clear delta theta alpha beta gammalow gammahigh broadband


pow_spec_corrmatrix = {};
pow_spec_corrmatrix.delta = delta_pow_corr_matrix;
pow_spec_corrmatrix.theta = theta_pow_corr_matrix;
pow_spec_corrmatrix.alpha = alpha_pow_corr_matrix;
pow_spec_corrmatrix.beta = beta_pow_corr_matrix;
pow_spec_corrmatrix.gammalow = gammalow_pow_corr_matrix;
pow_spec_corrmatrix.gammahigh = gammahigh_pow_corr_matrix;
pow_spec_corrmatrix.broadband = broadband_pow_corr_matrix;

clear delta_pow_corr_matrix theta_pow_corr_matrix
clear alpha_pow_corr_matrix beta_pow_corr_matrix
clear gammalow_pow_corr_matrix gammahigh_pow_corr_matrix
clear broadband_pow_corr_matrix

pow_spec_corr = {};
pow_spec_corr.delta = delta_pow_corr;
pow_spec_corr.theta = theta_pow_corr;
pow_spec_corr.alpha = alpha_pow_corr;
pow_spec_corr.beta = beta_pow_corr;
pow_spec_corr.gammalow = gammalow_pow_corr;
pow_spec_corr.gammahigh = gammahigh_pow_corr;
pow_spec_corr.broadband = broadband_pow_corr;

clear delta_pow_corr theta_pow_corr alpha_pow_corr beta_pow_corr
clear gammalow_pow_corr gammahigh_pow_corr broadband_pow_corr;

save(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/',subject,'/conn_analysis/',subject,'_run',run,'_icEEG_static_powspec_conn.mat'])