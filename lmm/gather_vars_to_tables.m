% script use to gather FC values into a table that canbe used to run LMM in R
for i = 1:length(subjects)
    sub = subjects(i,1);
    
    run = subjects(i,5);
    load(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/sub-0',num2str(sub),'/conn_analysis/sub-0',num2str(sub),'_run0',num2str(run),'_static_conn.mat'], 'static_conn')
    
    delta_stat = atanh(static_conn.delta_pow);
    theta_stat = atanh(static_conn.theta_pow);
    alpha_stat = atanh(static_conn.alpha_pow);
    beta_stat = atanh(static_conn.beta_pow);
    gammalow_stat = atanh(static_conn.gammalow_pow);
    gammahigh_stat = atanh(static_conn.gammahigh_pow);
    boldfc_stat = atanh(static_conn.bold_fc_static);
    
    
    id = [];
    id(1:length(delta_stat),1) = subjects(i,1); 

    sex = [];
    sex(1:length(delta_stat),1) = subjects(i,2); 

    age = [];
    age(1:length(delta_stat),1) = subjects(i,3); 

    mfd = [];
    mfd(1:length(delta_stat),1) = subjects(i,4); 

    contacts = [];
    contacts(1:length(delta_stat),1) = 1:length(delta_stat);
    
    data{i} = [id, sex, age, mfd, contacts, delta_stat, theta_stat, alpha_stat, beta_stat, gammalow_stat, gammahigh_stat, boldfc_stat];

% tbl = table(atanh(static_conn.delta_pow),atanh(static_conn.theta_pow),atanh(static_conn.alpha_pow),...
%     atanh(static_conn.beta_pow),atanh(static_conn.gammalow_pow),atanh(static_conn.gammahigh_pow),atanh(static_conn.bold_fc_static));

end

static = vertcat(data{1},data{2},data{3},data{4},data{5},data{6},data{7},data{8},...
    data{9},data{10},data{11},data{12},data{13},data{14},data{15},data{16},data{17});

static_data_longform = table(static(:,1),static(:,2),static(:,3),static(:,4),static(:,5),static(:,6),static(:,7),...
    static(:,8),static(:,9),static(:,10),static(:,11),static(:,12));

static_data_longform.Properties.VariableNames = {'ID','SEX','AGE','MFD','CONTACTS','DELTA_STAT','THETA_STAT','ALPHA_STAT','BETA_STAT','GAMMALOW_STAT','GAMMAHIGH_STAT','BOLDFC_STAT'};


%%
for i = 1:length(subjects)
    sub = subjects(i,1);
    
    run = subjects(i,5);
    load(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/sub-0',num2str(sub),'/conn_analysis/sliding_win/120s/sub-0',num2str(sub),'_run0',num2str(run),'_conn_metrics_120s.mat'],'conn_metrics_120s')
    load(['/Users/kristinasabaroedin/Documents/icEEG_fMRI/sub-0',num2str(sub),'/conn_analysis/sub-0',num2str(sub),'_run0',num2str(run),'_static_conn.mat'], 'static_conn')
    
 
    delta_stat = [];
    x = static_conn.delta_pow';
    delta_stat = repmat(x,1,18);
    delta_stat = atanh(delta_stat)';
    delta_stat = delta_stat(:);
    clear x
    
    theta_stat = [];
    x = static_conn.theta_pow';
    theta_stat = repmat(x,1,18);
    theta_stat = atanh(theta_stat)';
    theta_stat = theta_stat(:);
    clear x
    
    alpha_stat = [];
    x = static_conn.alpha_pow';
    alpha_stat = repmat(x,1,18);
    alpha_stat = atanh(alpha_stat)';
    alpha_stat = alpha_stat(:);
    clear x
    
    beta_stat = [];
    x = static_conn.beta_pow';
    beta_stat = repmat(x,1,18);
    beta_stat = atanh(beta_stat)';
    beta_stat = beta_stat(:);
    clear x
        
    gammalow_stat = [];
    x = static_conn.gammalow_pow';
    gammalow_stat = repmat(x,1,18);
    gammalow_stat = atanh(gammalow_stat)';
    gammalow_stat = gammalow_stat(:);
    clear x
    
    gammahigh_stat = [];
    x = static_conn.gammahigh_pow';
    gammahigh_stat = repmat(x,1,18);
    gammahigh_stat = atanh(gammahigh_stat)';
    gammahigh_stat = gammahigh_stat(:);
    clear x
            
    delta_dyn = atanh(conn_metrics_120s.powspec_corr_delta(:));
    theta_dyn = atanh(conn_metrics_120s.powspec_corr_theta(:));
    alpha_dyn = atanh(conn_metrics_120s.powspec_corr_alpha(:));
    beta_dyn = atanh(conn_metrics_120s.powspec_corr_beta(:));
    gammalow_dyn = atanh(conn_metrics_120s.powspec_corr_gammalow(:));
    gammahigh_dyn = atanh(conn_metrics_120s.powspec_corr_gammahigh(:));
    boldfc_dyn = atanh(conn_metrics_120s.boldfc(:));

    id = [];
    id(1:length(delta_dyn),1) = subjects(i,1); 

    sex = [];
    sex(1:length(delta_dyn),1) = subjects(i,2); 

    age = [];
    age(1:length(delta_dyn),1) = subjects(i,3); 

    mfd = [];
    mfd(1:length(delta_dyn),1) = subjects(i,4); 

    contacts = [];
    contacts(1:length(conn_metrics_120s.labelcmb),1) = 1:length(conn_metrics_120s.labelcmb);
    contacts = repmat(contacts,1,18);
    contacts = contacts(:);

    window = [];
    n=length(conn_metrics_120s.labelcmb); 
    x = (1:18)';
    window = repmat(x,1,n);
    window = window';
    window = window(:);
    
    delta_diff = delta_dyn - delta_stat;
    theta_diff = theta_dyn - theta_stat;
    alpha_diff = alpha_dyn - alpha_stat;
    beta_diff = beta_dyn - beta_stat;
    gammalow_diff = gammalow_dyn - gammalow_stat;
    gammahigh_diff = gammahigh_dyn - gammahigh_stat;
    
    delta_dyn_std = std(conn_metrics_120s.powspec_corr_delta')';
    delta_dyn_sd = [];
    delta_dyn_sd = repmat(delta_dyn_std,1,18);
    delta_dyn_sd = delta_dyn_sd(:);
    clear delta_dyn_std

    theta_dyn_std = std(conn_metrics_120s.powspec_corr_theta')';
    theta_dyn_sd = [];
    theta_dyn_sd = repmat(theta_dyn_std,1,18);
    theta_dyn_sd = theta_dyn_sd(:);
    clear theta_dyn_std

    alpha_dyn_std = std(conn_metrics_120s.powspec_corr_alpha')';
    alpha_dyn_sd = [];
    alpha_dyn_sd = repmat(alpha_dyn_std,1,18);
    alpha_dyn_sd = alpha_dyn_sd(:);
    clear alpha_dyn_std

    beta_dyn_std = std(conn_metrics_120s.powspec_corr_beta')';
    beta_dyn_sd = [];
    beta_dyn_sd = repmat(beta_dyn_std,1,18);
    beta_dyn_sd = beta_dyn_sd(:);
    clear beta_dyn_std

    gammalow_dyn_std = std(conn_metrics_120s.powspec_corr_gammalow')';
    gammalow_dyn_sd = [];
    gammalow_dyn_sd = repmat(gammalow_dyn_std,1,18);
    gammalow_dyn_sd = gammalow_dyn_sd(:);
    clear gammalow_dyn_std

    gammahigh_dyn_std = std(conn_metrics_120s.powspec_corr_gammahigh')';
    gammahigh_dyn_sd = [];
    gammahigh_dyn_sd = repmat(gammahigh_dyn_std,1,18);
    gammahigh_dyn_sd = gammahigh_dyn_sd(:);
    clear gammahigh_dyn_std

    data{i} = [id, sex, age, mfd, window, contacts, delta_stat, theta_stat, alpha_stat, beta_stat, gammalow_stat, gammahigh_stat,...
        delta_dyn, theta_dyn, alpha_dyn, beta_dyn, gammalow_dyn, gammahigh_dyn,...
        delta_diff, theta_diff, alpha_diff, beta_diff, gammalow_diff, gammahigh_diff, boldfc_dyn,...
        delta_dyn_sd, theta_dyn_sd, alpha_dyn_sd, beta_dyn_sd, gammalow_dyn_sd, gammahigh_dyn_sd];

end

dynamic = vertcat(data{1},data{2},data{3},data{4},data{5},data{6},data{7},data{8},...
    data{9},data{10},data{11},data{12},data{13},data{14},data{15},data{16},data{17});
dynamic_data_longform = table(dynamic(:,1),dynamic(:,2),dynamic(:,3),dynamic(:,4),dynamic(:,5),dynamic(:,6),dynamic(:,7),...
    dynamic(:,8),dynamic(:,9),dynamic(:,10),dynamic(:,11),dynamic(:,12),dynamic(:,13),dynamic(:,14),...
    dynamic(:,15),dynamic(:,16),dynamic(:,17),dynamic(:,18),dynamic(:,19),dynamic(:,20),dynamic(:,21),...
    dynamic(:,22),dynamic(:,23),dynamic(:,24),dynamic(:,25),dynamic(:,26),dynamic(:,27),dynamic(:,28),dynamic(:,29),dynamic(:,30),dynamic(:,31));
dynamic_data_longform.Properties.VariableNames = {'ID','SEX','AGE','MFD','WINDOW','CONTACTS','DELTA_STAT','THETA_STAT','ALPHA_STAT','BETA_STAT','GAMMALOW_STAT','GAMMAHIGH_STAT',...
    'DELTA_DYN','THETA_DYN','ALPHA_DYN','BETA_DYN','GAMMALOW_DYN','GAMMAHIGH_DYN',...
    'DELTA_DIFF','THETA_DIFF','ALPHA_DIFF','BETA_DIFF','GAMMALOW_DIFF','GAMMAHIGH_DIFF','BOLDFC_DYN',...
    'DELTA_DYN_SD','THETA_DYN_SD','ALPHA_DYN_SD','BETA_DYN_SD','GAMMALOW_DYN_SD','GAMMAHIGH_DYN_SD'};


filename = 'dyn_data_long_n17_02.xlsx';
writetable(dynamic_data_longform,filename,'Sheet',1)